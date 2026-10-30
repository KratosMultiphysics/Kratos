"""Process running a trained model inside a Kratos solution loop.

Per configured execution point, the process gathers the input fields into one
torch tensor of shape (n_entities, total_input_width), runs a no-grad forward
pass, splits the (n_entities, total_output_width) result by each output
field's width and writes the pieces back onto the model part.

Predictions are written into existing, physically-meaningful Kratos
variables, so downstream code cannot tell an ML-predicted value from a
solver-computed one. torch is imported lazily on first execution.
"""

import numpy

import KratosMultiphysics as Kratos
from KratosMultiphysics.PhysicsNeMoApplication.deployment import model_registry
from KratosMultiphysics.PhysicsNeMoApplication.deployment import ood_guard_utils
from KratosMultiphysics.PhysicsNeMoApplication.deployment import geometry_guard_utils
from KratosMultiphysics.PhysicsNeMoApplication.bridges import torch_bridge
from KratosMultiphysics.PhysicsNeMoApplication.utilities.tensor_adaptor_dataset_utils import GetTensorAdaptor
from KratosMultiphysics.PhysicsNeMoApplication.utilities.nvtx_utils import NvtxRange

_EXECUTION_POINTS = ("initialize_solution_step", "finalize_solution_step")


def GatherInputFields(model_part, field_specs, local_only: bool = False):
    """Gathers fields as a list of (n_entities, width) torch tensors.

    Shared by InferenceProcess and the CoSimulation surrogate wrapper.

    Args:
        local_only: If True, read owned entities only (the communicator's
            LocalMesh) rather than owned + ghost. Required on a distributed
            model part so the rows match CouplingInterfaceData's layout.

    Returns:
        (inputs, n_entities): the per-field tensors and their common entity
        count (all fields must live on the same entities).
    """
    inputs = []
    n_entities = None
    for variable_name, data_location in field_specs:
        variable = Kratos.KratosGlobals.GetVariable(variable_name)
        tensor_adaptor = GetTensorAdaptor(model_part, data_location, variable,
                                          local_only=local_only)
        data = torch_bridge.KratosTensorToTorch(tensor_adaptor)
        data = data.reshape(data.shape[0], -1)
        if n_entities is None:
            n_entities = data.shape[0]
        elif data.shape[0] != n_entities:
            raise ValueError(
                f"Input field \"{variable_name}\" has {data.shape[0]} entities but previous "
                f"inputs have {n_entities}; all input fields must live on the same entities.")
        inputs.append(data)
    return inputs, n_entities


def WriteOutputFields(model_part, field_specs, prediction, n_entities,
                      local_only: bool = False, normalization=None,
                      scale_only: bool = False) -> None:
    """Splits an (n_entities, total_width) prediction across the fields.

    With local_only the prediction is written to owned entities only; the
    caller is then responsible for synchronizing ghosts.

    normalization (a model_registry.LoadOutputNormalization entry) inverts
    the training normalization first, so a model trained on normalized
    targets writes physical values. None - the default - is the identity,
    which is what every configuration written before the card carried this
    does. Pass scale_only for a spread rather than a mean; see
    model_registry.ApplyOutputNormalization.
    """
    output_adaptors = []
    total_width = 0
    for variable_name, data_location in field_specs:
        variable = Kratos.KratosGlobals.GetVariable(variable_name)
        tensor_adaptor = GetTensorAdaptor(model_part, data_location, variable,
                                          local_only=local_only)
        width = int(numpy.prod(tensor_adaptor.data.shape[1:], dtype=int))
        output_adaptors.append((tensor_adaptor, width))
        total_width += width

    if normalization is not None:
        # before the shape check, so a card/model mismatch is reported
        # against the widths the fields actually require
        prediction = model_registry.ApplyOutputNormalization(
            prediction, normalization, scale_only=scale_only)

    if tuple(prediction.shape) != (n_entities, total_width):
        raise ValueError(
            f"Model returned shape {list(prediction.shape)} but the configured output "
            f"fields require ({n_entities}, {total_width}).")

    offset = 0
    for tensor_adaptor, width in output_adaptors:
        chunk = prediction[:, offset:offset + width].reshape(tensor_adaptor.data.shape)
        torch_bridge.TorchToKratosTensor(chunk, tensor_adaptor)
        offset += width


def SynchronizeOutputFields(model_part, field_specs) -> None:
    """Refreshes every rank's ghost copies after an owned-only write.

    The counterpart of WriteOutputFields(local_only=True): the owned values
    are authoritative, so each ghost is overwritten from its owner. A no-op
    on a serial model part. Shared by the CoSimulation surrogate wrapper and
    the domain-parallel write-back.
    """
    if not model_part.IsDistributed():
        return
    communicator = model_part.GetCommunicator()
    for variable_name, data_location in field_specs:
        variable = Kratos.KratosGlobals.GetVariable(variable_name)
        if hasattr(variable, "GetSourceVariable"):
            variable = variable.GetSourceVariable()
        if data_location == "node_historical":
            communicator.SynchronizeVariable(variable)
        elif data_location == "node_non_historical":
            communicator.SynchronizeNonHistoricalVariable(variable)


def Factory(settings: Kratos.Parameters, model: Kratos.Model) -> "InferenceProcess":
    if not isinstance(settings, Kratos.Parameters):
        raise TypeError("Expected input shall be a Parameters object, encapsulating a json string")
    return InferenceProcess(model, settings["Parameters"])


class InferenceProcess(Kratos.Process):
    """Runs model inference each output_interval steps at an execution point."""

    def __init__(self, model: Kratos.Model, settings: Kratos.Parameters) -> None:
        super().__init__()

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
                    "data_location" : "node_non_historical"
                }
            ],
            "execution_point" : "finalize_solution_step",
            "output_interval" : 1,
            "ood_guard"       : {},
            "geometry_guard"  : {},
            "uncertainty"     : {}
        }""")
        settings.ValidateAndAssignDefaults(default_settings)
        for key in ("input_fields", "output_fields"):
            for i in range(settings[key].size()):
                settings[key][i].ValidateAndAssignDefaults(default_settings[key][0])

        # optional OOD guard on the gathered inputs (see ood_guard_utils)
        self._ood_guard = ood_guard_utils.GuardCheck(settings["ood_guard"])
        # and on the SHAPE, which the field guard cannot see at all
        self._geometry_guard = geometry_guard_utils.GeometryGuardCheck(
            settings["geometry_guard"])

        # optional predictive uncertainty (MC dropout or checkpoint ensemble)
        uncertainty_defaults = Kratos.Parameters("""{
            "method"             : "none",
            "num_samples"        : 16,
            "seed"               : -1,
            "gp_head_file"       : "",
            "retain_ensemble"    : false,
            "gp_feature_fields"  : [],
            "uncertainty_fields" : []
        }""")
        uncertainty = settings["uncertainty"]
        uncertainty.ValidateAndAssignDefaults(uncertainty_defaults)
        for i in range(uncertainty["uncertainty_fields"].size()):
            uncertainty["uncertainty_fields"][i].ValidateAndAssignDefaults(
                default_settings["output_fields"][0])
        self.uncertainty_method = uncertainty["method"].GetString()
        self.uncertainty_samples = uncertainty["num_samples"].GetInt()
        self.uncertainty_seed = uncertainty["seed"].GetInt()
        self.uncertainty_specs = self._ReadFieldSpecs(uncertainty["uncertainty_fields"])
        self.gp_head_file = uncertainty["gp_head_file"].GetString()
        # The (M, ...) stack below is built and then reduced to mean/std; the
        # members themselves are the only thing a proper scoring rule like
        # CRPS can consume, so they are optionally kept for a validation
        # process to read (the ValidationMetricsProcess.history pattern).
        self.retain_ensemble = uncertainty["retain_ensemble"].GetBool()
        self.last_ensemble = None
        for i in range(uncertainty["gp_feature_fields"].size()):
            uncertainty["gp_feature_fields"][i].ValidateAndAssignDefaults(
                default_settings["input_fields"][0])
        self.gp_feature_specs = self._ReadFieldSpecs(uncertainty["gp_feature_fields"])
        self._gp_head = None
        if self.uncertainty_method not in ("none", "mc_dropout", "ensemble", "gp"):
            raise ValueError(
                f"Unsupported uncertainty method \"{self.uncertainty_method}\". "
                "Use \"none\", \"mc_dropout\", \"ensemble\" or \"gp\".")
        if self.uncertainty_method == "gp" and not self.gp_head_file:
            raise ValueError(
                "The \"gp\" uncertainty method needs \"gp_head_file\": the fitted head "
                "is a sidecar next to the checkpoint (see uncertainty_utils.SaveGpHead).")

        self.model_part = model[settings["model_part_name"].GetString()]
        self.model_settings = settings["model_settings"].Clone()
        self.input_specs = self._ReadFieldSpecs(settings["input_fields"])
        self.output_specs = self._ReadFieldSpecs(settings["output_fields"])
        self.execution_point = settings["execution_point"].GetString()
        self.output_interval = settings["output_interval"].GetInt()

        if self.execution_point not in _EXECUTION_POINTS:
            raise ValueError(
                f"Unsupported execution point \"{self.execution_point}\". "
                f"Supported: {', '.join(_EXECUTION_POINTS)}.")
        if self.output_interval < 1:
            raise ValueError(f"\"output_interval\" must be >= 1 [ output_interval = {self.output_interval} ].")

        # Model loading is deferred to the first execution so that merely
        # constructing the process never requires torch.
        self._model = None
        self._device = None
        self._ensemble = None

    @staticmethod
    def _ReadFieldSpecs(fields: Kratos.Parameters):
        return [
            (fields[i]["variable_name"].GetString(), fields[i]["data_location"].GetString())
            for i in range(fields.size())
        ]

    def ExecuteInitializeSolutionStep(self) -> None:
        if self.execution_point == "initialize_solution_step":
            self._ExecuteIfDue()

    def ExecuteFinalizeSolutionStep(self) -> None:
        if self.execution_point == "finalize_solution_step":
            self._ExecuteIfDue()

    def _ExecuteIfDue(self) -> None:
        if self.model_part.ProcessInfo[Kratos.STEP] % self.output_interval == 0:
            self.RunInference()

    def _GetModel(self):
        if self._model is None:
            self._model, self._device = model_registry.LoadModelWithCardCheck(
                self.model_settings, self.input_specs, self.output_specs, type(self).__name__)
        return self._model

    def _NormalizationCardFile(self):
        """The checkpoint path the model card is keyed off, or None to let
        the registry resolve "checkpoint_file"/"checkpoint_files" itself.
        ONNX and Triton override this: they key off "onnx_file" and
        "card_file", where the base lookup would find nothing."""
        return None

    def _GetNormalization(self):
        """The card's output de-normalization, or None.

        Loaded lazily HERE rather than as a side effect of _GetModel: an
        ensemble deployment never calls _GetModel at all (it goes through
        _GetEnsembleModels), and hanging the load off one of the two paths
        meant the other silently skipped de-normalization entirely.
        """
        if not hasattr(self, "_normalization"):
            self._normalization = model_registry.LoadOutputNormalization(
                self.model_settings, checkpoint_file=self._NormalizationCardFile())
        return self._normalization

    def _GetInputNormalization(self):
        """The card's input normalization, or None - the symmetric half of
        _GetNormalization (see model_registry.LoadInputNormalization)."""
        if not hasattr(self, "_input_normalization"):
            self._input_normalization = model_registry.LoadInputNormalization(
                self.model_settings, checkpoint_file=self._NormalizationCardFile())
        return self._input_normalization

    def _GatherFeatures(self):
        """The input fields as ONE (n_entities, total_width) tensor, with the
        card's "input_normalization" applied.

        Every subclass feeds its model from here, so the key reaches all of
        them at once - and before the OOD guard check, which was calibrated
        on what the model saw in training.

        Returns:
            (features, n_entities).
        """
        torch = torch_bridge._TryImportTorch()
        inputs, n_entities = self._GatherInputs()
        features = torch.cat(inputs, dim=-1)
        features = model_registry.ApplyInputNormalization(features, self._GetInputNormalization())
        return features, n_entities

    def _GetEnsembleModels(self):
        """Loads the checkpoint ensemble named by model_settings["checkpoint_files"]."""
        if self._ensemble is None:
            if not self.model_settings.Has("checkpoint_files"):
                raise ValueError(
                    "Ensemble uncertainty needs a \"checkpoint_files\" string array in "
                    "\"model_settings\" naming the member checkpoints.")
            files = self.model_settings["checkpoint_files"].GetStringArray()
            if len(files) < 2:
                raise ValueError(
                    f"An ensemble needs at least 2 checkpoints; got {len(files)}.")
            member_settings = self.model_settings.Clone()
            member_settings.RemoveValue("checkpoint_files")
            if not member_settings.Has("checkpoint_file"):
                member_settings.AddEmptyValue("checkpoint_file").SetString("")
            models = []
            for checkpoint_file in files:
                settings = member_settings.Clone()
                settings["checkpoint_file"].SetString(checkpoint_file)
                model, self._device = model_registry.LoadModelWithCardCheck(
                    settings, self.input_specs, self.output_specs, type(self).__name__)
                models.append(model)
            self._ensemble = models
        return self._ensemble

    def _CheckOOD(self, features) -> None:
        self._ood_guard.Check(features, type(self).__name__)
        self._CheckGeometry()

    def _CheckGeometry(self) -> None:
        """The shape guardrail, re-evaluated only when the mesh changes."""
        self._geometry_guard.Check(self.model_part, type(self).__name__)

    def _PredictWithUncertainty(self, forward_fn):
        """Runs forward_fn per the configured uncertainty method.

        forward_fn: callable model -> (n_entities, total_width) tensor.

        Returns:
            (prediction, std): the (mean) prediction and the per-entity
            standard deviation (None for method "none").
        """
        if self.uncertainty_method == "mc_dropout":
            from KratosMultiphysics.PhysicsNeMoApplication.deployment import uncertainty_utils
            return uncertainty_utils.MonteCarloPredict(
                self._GetModel(), forward_fn, self.uncertainty_samples, self.uncertainty_seed)
        if self.uncertainty_method == "ensemble":
            torch = torch_bridge._TryImportTorch()
            samples = torch.stack([forward_fn(model) for model in self._GetEnsembleModels()])
            if self.retain_ensemble:
                self.last_ensemble = samples.detach().clone()
            return samples.mean(dim=0), samples.std(dim=0)
        if self.uncertainty_method == "gp":
            return self._PredictWithGpHead(forward_fn)
        return forward_fn(self._GetModel()), None

    def _PredictWithGpHead(self, forward_fn):
        """Prediction plus the GP head's calibrated standard deviation.

        The head is fitted on the backbone's features, so it needs those
        features - "gp_feature_fields" names them, defaulting to the
        process' own input fields when it is left empty.
        """
        from KratosMultiphysics.PhysicsNeMoApplication.deployment import uncertainty_utils
        torch = torch_bridge._TryImportTorch()
        prediction = forward_fn(self._GetModel())

        if self.gp_feature_specs:
            features, _ = GatherInputFields(self.model_part, self.gp_feature_specs)
            features = torch.cat(features, dim=-1)
        else:  # the model's own inputs, normalized exactly as the model sees them
            features, _ = self._GatherFeatures()

        if self._gp_head is None:
            self._gp_head = uncertainty_utils.LoadGpHead(
                self.gp_head_file, input_dim=int(features.shape[-1]))
        _, std = uncertainty_utils.PredictWithGpHead(self._gp_head, features)

        if std.shape[-1] == 1 and prediction.shape[-1] > 1:
            std = std.expand_as(prediction)
        return prediction, std.to(prediction.dtype)

    def _WriteUncertainty(self, std, n_entities) -> None:
        if std is None or not self.uncertainty_specs:
            return
        # scale_only: this is a SPREAD. De-normalizing it like a mean would
        # shift every standard deviation by the training mean, which is
        # meaningless - the single easiest mistake for a shared hook.
        WriteOutputFields(self.model_part, self.uncertainty_specs, std, n_entities,
                          normalization=self._GetNormalization(), scale_only=True)

    def _GatherInputs(self):
        """Gathers the input fields as a list of (n_entities, width) tensors."""
        return GatherInputFields(self.model_part, self.input_specs)

    def _WriteOutputs(self, prediction, n_entities) -> None:
        """Splits an (n_entities, total_width) prediction across the output fields."""
        WriteOutputFields(self.model_part, self.output_specs, prediction, n_entities,
                          normalization=self._GetNormalization())

    def RunInference(self) -> None:
        torch = torch_bridge._TryImportTorch()

        with NvtxRange("PhysicsNeMo::GatherInputs"):
            features, n_entities = self._GatherFeatures()
        self._CheckOOD(features)

        def forward(model):
            with torch.no_grad(), NvtxRange("PhysicsNeMo::Forward"):
                return model(features.to(self._device)).cpu()

        prediction, std = self._PredictWithUncertainty(forward)
        with NvtxRange("PhysicsNeMo::WriteOutputs"):
            self._WriteOutputs(prediction, n_entities)
            self._WriteUncertainty(std, n_entities)
