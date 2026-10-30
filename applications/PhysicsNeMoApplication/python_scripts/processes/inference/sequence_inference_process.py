"""Process deploying one-to-many grid sequence models (RNN surrogates).

The deployment counterpart of physicsnemo's RNN pattern (One2ManyRNN /
Seq2SeqRNN, the 2D Navier-Stokes / Gray-Scott examples): at the first due
execution the configured input fields are sampled onto a regular grid and
fed to the model ONCE as the initial state - the model returns the whole
predicted future, (1, C, 1, *spatial) -> (1, C, T, *spatial). Each
subsequent due execution pops the next predicted state from the buffer and
scatters it onto the output fields, until the T predicted steps are
exhausted (further executions warn once and do nothing).

Planar 2D cases use the thin-axis idiom (see grid_dataset_export_process):
"squeeze_axis" collapses the thin spatial axis by its mean before the
forward pass (dimension=2 models) and duplicates the prediction across it
on the way back.

"window_as_time_axis" targets spatiotemporal block operators - FNO with
dimension=4 (its time axis is just a fourth grid axis) and seq2seq RNNs:
instead of seeding from a single state, the process accumulates the sampled
grid at each due step into a rolling window of "window_size" states
(warm-up steps are logged, nothing written) and, once full, feeds the model
the whole (1, C, K, *spatial) block; the returned (C, T, *spatial) block is
buffered as usual. FNO(dimension=4) returns T == K - the predicted NEXT
block of K states, the standard FNO time-block surrogate.

"model_interface" : "dpot" deploys physicsnemo.models.dpot.DPOTNet, the
AFNO-mixing PDE foundation model. It is an ordinary grid series surrogate
here - it reads a window of states and returns the next ones - but its
tensor layout is its own, (B, H, W, T, C) rather than the (B, C, T, H, W)
every other grid model in this application uses, so the interface is a pair
of permutations. Pretrained DPOT weights are not verified from here; the
architecture and the fine-tuning path (domino_finetune's LoRA helpers apply
unchanged, with a pattern matching DPOT's blocks) are.

torch is imported lazily on first execution.
"""

import numpy

import KratosMultiphysics as Kratos
from KratosMultiphysics.PhysicsNeMoApplication.bridges import grid_bridge
from KratosMultiphysics.PhysicsNeMoApplication.deployment import model_registry
from KratosMultiphysics.PhysicsNeMoApplication.bridges import torch_bridge
_EXECUTION_POINTS = ("initialize_solution_step", "finalize_solution_step")
_MODEL_INTERFACES = ("grid", "dpot")


def Factory(settings: Kratos.Parameters, model: Kratos.Model) -> "SequenceInferenceProcess":
    if not isinstance(settings, Kratos.Parameters):
        raise TypeError("Expected input shall be a Parameters object, encapsulating a json string")
    return SequenceInferenceProcess(model, settings["Parameters"])


class SequenceInferenceProcess(Kratos.Process):
    """Seeds a sequence model once, then writes one predicted state per due step."""

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
                    "data_location" : "node_historical"
                }
            ],
            "grid_shape"          : [8, 8, 8],
            "bounding_box"        : [],
            "squeeze_axis"        : -1,
            "model_interface"     : "grid",
            "window_as_time_axis" : false,
            "window_size"         : 2,
            "execution_point"     : "finalize_solution_step",
            "output_interval"     : 1
        }""")
        settings.ValidateAndAssignDefaults(default_settings)
        for key in ("input_fields", "output_fields"):
            for i in range(settings[key].size()):
                settings[key][i].ValidateAndAssignDefaults(default_settings[key][0])

        self.model_part = model[settings["model_part_name"].GetString()]
        self.model_settings = settings["model_settings"].Clone()
        self.input_specs = self._ReadFieldSpecs(settings["input_fields"])
        self.output_specs = self._ReadFieldSpecs(settings["output_fields"])
        self.grid_shape = tuple(int(n) for n in settings["grid_shape"].GetVector())
        if len(self.grid_shape) != 3:
            raise ValueError(f"\"grid_shape\" must have three entries, got {self.grid_shape}.")

        box = settings["bounding_box"].GetVector()
        if len(box) == 0:
            self.bounding_box = None  # resolved at the seeding execution
        elif len(box) == 6:
            self.bounding_box = (numpy.array(box[:3]), numpy.array(box[3:]))
        else:
            raise ValueError("\"bounding_box\" must be empty or [x0,y0,z0,x1,y1,z1].")

        squeeze_axis = settings["squeeze_axis"].GetInt()
        if squeeze_axis == -1:
            self.squeeze_axis = None
        elif squeeze_axis in (0, 1, 2):
            self.squeeze_axis = squeeze_axis
        else:
            raise ValueError(f"\"squeeze_axis\" must be -1 (off), 0, 1 or 2, got {squeeze_axis}.")

        self.model_interface = settings["model_interface"].GetString()
        if self.model_interface not in _MODEL_INTERFACES:
            raise ValueError(
                f"Unsupported model interface \"{self.model_interface}\". Supported: "
                f"{', '.join(_MODEL_INTERFACES)}.")
        self.window_as_time_axis = settings["window_as_time_axis"].GetBool()
        self.window_size = settings["window_size"].GetInt()
        if self.window_as_time_axis and self.window_size < 2:
            raise ValueError(f"\"window_size\" must be >= 2 [ window_size = {self.window_size} ].")

        self.execution_point = settings["execution_point"].GetString()
        if self.execution_point not in _EXECUTION_POINTS:
            raise ValueError(
                f"Unsupported execution point \"{self.execution_point}\". "
                f"Supported: {', '.join(_EXECUTION_POINTS)}.")
        self.output_interval = settings["output_interval"].GetInt()
        if self.output_interval < 1:
            raise ValueError(f"\"output_interval\" must be >= 1 [ output_interval = {self.output_interval} ].")

        self._model = None
        self._normalization = None
        self._device = None
        self._buffer = None  # None = not seeded yet; list afterwards
        self._window = []    # window_as_time_axis: sampled states awaiting seeding
        self._exhaustion_warned = False

    @staticmethod
    def _ReadFieldSpecs(fields: Kratos.Parameters):
        return [
            (fields[i]["variable_name"].GetString(), fields[i]["data_location"].GetString())
            for i in range(fields.size())
        ]

    @property
    def predicted_steps_left(self) -> int:
        return len(self._buffer) if self._buffer is not None else 0

    def ExecuteInitializeSolutionStep(self) -> None:
        if self.execution_point == "initialize_solution_step":
            self._ExecuteIfDue()

    def ExecuteFinalizeSolutionStep(self) -> None:
        if self.execution_point == "finalize_solution_step":
            self._ExecuteIfDue()

    def _ExecuteIfDue(self) -> None:
        if self.model_part.ProcessInfo[Kratos.STEP] % self.output_interval != 0:
            return
        if self._buffer is None:
            if self.window_as_time_axis:
                self._window.append(self._SampleGrid())
                if len(self._window) < self.window_size:
                    Kratos.Logger.PrintInfo(
                        type(self).__name__,
                        f"Accumulating the input window ({len(self._window)}/"
                        f"{self.window_size}); no prediction yet.")
                    return
            self._SeedRollout()
        else:
            self._WriteNextState()

    def _GetInputNormalization(self):
        """The card's input normalization, or None. Lazy rather than loaded
        with the model: the input window accumulates BEFORE the model is
        loaded, and every sampled grid must already be normalized."""
        if not hasattr(self, "_input_normalization"):
            self._input_normalization = model_registry.LoadInputNormalization(self.model_settings)
        return self._input_normalization

    def _SampleGrid(self):
        if self.bounding_box is None:
            self.bounding_box = grid_bridge.ComputeBoundingBox(self.model_part)
        grid, _ = grid_bridge.SampleFieldsOnGrid(
            self.model_part, self.input_specs, self.grid_shape, self.bounding_box)
        if self.squeeze_axis is not None:
            grid = grid.mean(axis=1 + self.squeeze_axis)
        # the card's "input_normalization", per channel along axis 0
        return model_registry.ApplyInputNormalization(
            grid, self._GetInputNormalization(), channel_axis=0)

    def _SeedRollout(self) -> None:
        torch = torch_bridge._TryImportTorch()
        if self._model is None:
            self._model, self._device = model_registry.LoadModelWithCardCheck(
                self.model_settings, self.input_specs, self.output_specs, type(self).__name__)
            self._normalization = model_registry.LoadOutputNormalization(self.model_settings)

        if self.window_as_time_axis:
            # (1, C, K, *spatial): the accumulated window as the time axis
            grid = self._window[0]
            batch_array = numpy.stack(self._window, axis=1)[None]
            self._window = []
        else:
            grid = self._SampleGrid()
            batch_array = grid[None, :, None]  # (1, C, 1, *spatial)

        with torch.no_grad():
            batch = torch.from_numpy(batch_array).to(self._device)
            parameter = next(self._model.parameters(), None)
            if parameter is not None:
                batch = batch.to(parameter.dtype)
            if self.model_interface == "dpot":
                # DPOT reads (B, H, W, T, C) - the time axis LAST but one and
                # the channels last - where every other grid model here reads
                # (B, C, T, H, W). The permutation is the whole adapter.
                if batch.ndim != 5:
                    raise ValueError(
                        "The \"dpot\" interface consumes a 2-D grid series shaped "
                        f"(1, C, T, H, W); got {tuple(batch.shape)}. Planar cases reach "
                        "that through \"squeeze_axis\", and \"window_as_time_axis\" "
                        "supplies the time axis.")
                expected = getattr(self._model, "in_timesteps", None)
                if expected is not None and int(expected) != int(batch.shape[2]):
                    raise ValueError(
                        f"The model was built for in_timesteps = {int(expected)} but the "
                        f"window supplies {int(batch.shape[2])}; set \"window_size\" to "
                        "match.")
                prediction = self._model(
                    batch.permute(0, 3, 4, 2, 1)).permute(0, 4, 3, 1, 2)
                prediction = prediction.cpu().numpy()[0]  # (C_out, T_out, H, W)
            else:
                prediction = self._model(batch).cpu().numpy()[0]  # (C, T, *spatial)

        if prediction.ndim != grid.ndim + 1:
            raise ValueError(
                f"The model must return a (C, T, *spatial) sequence per sample; got shape "
                f"{prediction.shape} for a {grid.ndim - 1}-dimensional spatial grid.")
        self._buffer = [prediction[:, t] for t in range(prediction.shape[1])]
        Kratos.Logger.PrintInfo(
            type(self).__name__,
            f"Seeded the rollout at step {self.model_part.ProcessInfo[Kratos.STEP]}: "
            f"{len(self._buffer)} predicted state(s) buffered.")

    def _WriteNextState(self) -> None:
        if not self._buffer:
            if not self._exhaustion_warned:
                self._exhaustion_warned = True
                Kratos.Logger.PrintWarning(
                    type(self).__name__,
                    "All predicted states have been consumed; further steps write nothing. "
                    "Increase the model's nr_tsteps (or re-seed with a new process) for "
                    "longer rollouts.")
            return

        state = numpy.asarray(self._buffer.pop(0), dtype=numpy.float64)
        # the card's "output_normalization" makes a normalized state physical
        # (per channel along axis 0, the grid layout)
        state = model_registry.ApplyOutputNormalization(state, self._normalization, channel_axis=0)
        if self.squeeze_axis is not None:
            # duplicate the prediction across the collapsed thin axis
            thin_size = self.grid_shape[self.squeeze_axis]
            state = numpy.repeat(
                numpy.expand_dims(state, 1 + self.squeeze_axis), thin_size, axis=1 + self.squeeze_axis)
        grid_bridge.ScatterGridToNodes(state, self.bounding_box, self.model_part, self.output_specs)
