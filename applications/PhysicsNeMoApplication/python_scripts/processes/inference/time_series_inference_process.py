"""Process running an autoregressive time-series surrogate in a solution loop.

For transient problems: the process keeps a rolling history of the last K
gathered input states and, once the history is full, predicts the next state
node-locally each step. Model contract: input (N, K*W_in) — the history
concatenated along channels, OLDEST FIRST — output (N, W_out).

The history is appended at every execution point call (before predicting),
so during the first K-1 steps the process only warms up (logged, no write).

torch is imported lazily on first prediction.
"""

import collections

import numpy

import KratosMultiphysics as Kratos
from KratosMultiphysics.PhysicsNeMoApplication.deployment import model_registry
from KratosMultiphysics.PhysicsNeMoApplication.deployment import ood_guard_utils
from KratosMultiphysics.PhysicsNeMoApplication.bridges import torch_bridge
from KratosMultiphysics.PhysicsNeMoApplication.utilities.tensor_adaptor_dataset_utils import GetTensorAdaptor

_EXECUTION_POINTS = ("initialize_solution_step", "finalize_solution_step")


def Factory(settings: Kratos.Parameters, model: Kratos.Model) -> "TimeSeriesInferenceProcess":
    if not isinstance(settings, Kratos.Parameters):
        raise TypeError("Expected input shall be a Parameters object, encapsulating a json string")
    return TimeSeriesInferenceProcess(model, settings["Parameters"])


class TimeSeriesInferenceProcess(Kratos.Process):
    """Autoregressive next-state prediction from a K-state history."""

    def __init__(self, model: Kratos.Model, settings: Kratos.Parameters) -> None:
        super().__init__()

        # optional OOD guard on the input window (see ood_guard_utils)
        ood_settings = Kratos.Parameters("{}")
        if settings.Has("ood_guard"):
            ood_settings = settings["ood_guard"].Clone()
            settings.RemoveValue("ood_guard")
        self._ood_guard = ood_guard_utils.GuardCheck(ood_settings)

        default_settings = Kratos.Parameters("""{
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "model_settings"  : {},
            "input_fields"    : [
                {
                    "variable_name" : "PLEASE_SPECIFY_VARIABLE_NAME",
                    "data_location" : "node_historical"
                }
            ],
            "output_fields"   : [
                {
                    "variable_name" : "PLEASE_SPECIFY_VARIABLE_NAME",
                    "data_location" : "node_historical"
                }
            ],
            "history_size"    : 2,
            "execution_point" : "finalize_solution_step",
            "output_interval" : 1
        }""")
        settings.ValidateAndAssignDefaults(default_settings)
        for key in ("input_fields", "output_fields"):
            for i in range(settings[key].size()):
                settings[key][i].ValidateAndAssignDefaults(default_settings[key][0])

        self.model_part = model[settings["model_part_name"].GetString()]
        self.model_settings = settings["model_settings"].Clone()
        self.input_specs = [
            (settings["input_fields"][i]["variable_name"].GetString(),
             settings["input_fields"][i]["data_location"].GetString())
            for i in range(settings["input_fields"].size())
        ]
        self.output_specs = [
            (settings["output_fields"][i]["variable_name"].GetString(),
             settings["output_fields"][i]["data_location"].GetString())
            for i in range(settings["output_fields"].size())
        ]
        self.history_size = settings["history_size"].GetInt()
        if self.history_size < 2:
            raise ValueError(f"\"history_size\" must be >= 2 [ history_size = {self.history_size} ].")
        self.execution_point = settings["execution_point"].GetString()
        if self.execution_point not in _EXECUTION_POINTS:
            raise ValueError(
                f"Unsupported execution point \"{self.execution_point}\". "
                f"Supported: {', '.join(_EXECUTION_POINTS)}.")
        self.output_interval = settings["output_interval"].GetInt()
        if self.output_interval < 1:
            raise ValueError(f"\"output_interval\" must be >= 1 [ output_interval = {self.output_interval} ].")

        self.history = collections.deque(maxlen=self.history_size)
        self._model = None
        self._device = None
        self._normalization = None
        self._input_normalization = None

    def ExecuteInitializeSolutionStep(self) -> None:
        if self.execution_point == "initialize_solution_step":
            self._ExecuteIfDue()

    def ExecuteFinalizeSolutionStep(self) -> None:
        if self.execution_point == "finalize_solution_step":
            self._ExecuteIfDue()

    def _GatherInputs(self) -> numpy.ndarray:
        fields = []
        for variable_name, data_location in self.input_specs:
            variable = Kratos.KratosGlobals.GetVariable(variable_name)
            tensor_adaptor = GetTensorAdaptor(self.model_part, data_location, variable)
            data = numpy.array(tensor_adaptor.data)
            fields.append(data.reshape(data.shape[0], -1))
        return numpy.concatenate(fields, axis=1)

    def _ExecuteIfDue(self) -> None:
        if self.model_part.ProcessInfo[Kratos.STEP] % self.output_interval != 0:
            return
        self.history.append(self._GatherInputs())
        if len(self.history) < self.history_size:
            Kratos.Logger.PrintInfo(
                "TimeSeriesInferenceProcess",
                f"Warming up history ({len(self.history)}/{self.history_size}); no prediction yet.")
            return
        self.RunInference()

    def RunInference(self) -> None:
        torch = torch_bridge._TryImportTorch()
        if self._model is None:
            self._model, self._device = model_registry.LoadModelWithCardCheck(
                self.model_settings, self.input_specs, self.output_specs, type(self).__name__)
            self._normalization = model_registry.LoadOutputNormalization(self.model_settings)
            self._input_normalization = model_registry.LoadInputNormalization(self.model_settings)

        window = numpy.concatenate(list(self.history), axis=1)  # (N, K*W_in), oldest first
        # the card's "input_normalization" (width K*W_in, the whole window)
        window = model_registry.ApplyInputNormalization(window, self._input_normalization)
        if self._ood_guard.enabled:
            self._ood_guard.Check(torch.from_numpy(window), type(self).__name__)
        parameter = next(self._model.parameters(), None)
        dtype = parameter.dtype if parameter is not None else torch.float64
        with torch.no_grad():
            prediction = self._model(
                torch.from_numpy(window).to(self._device, dtype)).cpu().double().numpy()
        # the card's "output_normalization" makes a normalized prediction physical
        prediction = model_registry.ApplyOutputNormalization(prediction, self._normalization)

        offset = 0
        for variable_name, data_location in self.output_specs:
            variable = Kratos.KratosGlobals.GetVariable(variable_name)
            tensor_adaptor = GetTensorAdaptor(self.model_part, data_location, variable)
            width = int(numpy.prod(tensor_adaptor.data.shape[1:], dtype=int))
            tensor_adaptor.data[:] = prediction[:, offset:offset + width].reshape(tensor_adaptor.data.shape)
            tensor_adaptor.StoreData()
            offset += width
        if offset != prediction.shape[1]:
            raise ValueError(
                f"Model returned {prediction.shape[1]} channels but the output fields consume {offset}.")
