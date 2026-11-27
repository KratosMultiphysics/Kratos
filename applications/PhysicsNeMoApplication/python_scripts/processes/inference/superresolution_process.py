"""Process superresolving coarse-mesh fields onto a fine mesh.

The bridge for physicsnemo's voxel-grid superresolution models (e.g.
physicsnemo.models.srrn.SRResNet): per execution, the configured input
fields of the coarse model part are sampled onto a regular (C, D, H, W)
grid, the model upscales it by its scaling factor, and the fine grid is
scattered (trilinearly) onto the fine model part's nodes.

Planar 2D cases use the thin-axis idiom (see grid_dataset_export_process):
"squeeze_axis" collapses the thin spatial axis by its mean before the
forward pass - the model then sees (B, C, A, B') 2D grids, the layout of
physicsnemo's 2D operators (FNO dimension=2, AFNO/ModAFNO, DLWP, 2D
UNets) - and the prediction is duplicated across the thin axis on the way
back. Time-modulated models deploy through "model_interface": "modafno",
which passes the model part's TIME as the second (timestep) input,
matching ModAFNO.forward(x, mod).

torch/physicsnemo are imported lazily on first execution.
"""

import numpy

import KratosMultiphysics as Kratos
from KratosMultiphysics.PhysicsNeMoApplication.bridges import grid_bridge
from KratosMultiphysics.PhysicsNeMoApplication.deployment import model_registry
from KratosMultiphysics.PhysicsNeMoApplication.deployment import ood_guard_utils
from KratosMultiphysics.PhysicsNeMoApplication.bridges import torch_bridge
from KratosMultiphysics.PhysicsNeMoApplication.utilities.nvtx_utils import NvtxRange

_EXECUTION_POINTS = ("initialize_solution_step", "finalize_solution_step")


def Factory(settings: Kratos.Parameters, model: Kratos.Model) -> "SuperResolutionProcess":
    if not isinstance(settings, Kratos.Parameters):
        raise TypeError("Expected input shall be a Parameters object, encapsulating a json string")
    return SuperResolutionProcess(model, settings["Parameters"])


class SuperResolutionProcess(Kratos.Process):
    """Runs grid superresolution every output_interval steps."""

    def __init__(self, model: Kratos.Model, settings: Kratos.Parameters) -> None:
        super().__init__()

        # optional OOD guard on the sampled input grid (see ood_guard_utils)
        ood_settings = Kratos.Parameters("{}")
        if settings.Has("ood_guard"):
            ood_settings = settings["ood_guard"].Clone()
            settings.RemoveValue("ood_guard")
        self._ood_guard = ood_guard_utils.GuardCheck(ood_settings)

        default_settings = Kratos.Parameters("""{
            "coarse_model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "fine_model_part_name"   : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "model_settings"         : {},
            "input_fields"           : [
                {
                    "variable_name" : "PLEASE_SPECIFY_VARIABLE_NAME",
                    "data_location" : "node_historical"
                }
            ],
            "output_fields"          : [
                {
                    "variable_name" : "PLEASE_SPECIFY_VARIABLE_NAME",
                    "data_location" : "node_historical"
                }
            ],
            "coarse_grid_shape"      : [8, 8, 8],
            "bounding_box"           : [],
            "squeeze_axis"           : -1,
            "model_interface"        : "grid",
            "execution_point"        : "finalize_solution_step",
            "output_interval"        : 1
        }""")
        settings.ValidateAndAssignDefaults(default_settings)
        for key in ("input_fields", "output_fields"):
            for i in range(settings[key].size()):
                settings[key][i].ValidateAndAssignDefaults(default_settings[key][0])

        self.coarse_model_part = model[settings["coarse_model_part_name"].GetString()]
        self.fine_model_part = model[settings["fine_model_part_name"].GetString()]
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
        self.coarse_grid_shape = tuple(settings["coarse_grid_shape"].GetVector())
        self.coarse_grid_shape = tuple(int(n) for n in self.coarse_grid_shape)
        if len(self.coarse_grid_shape) != 3:
            raise ValueError(f"\"coarse_grid_shape\" must have three entries, got {self.coarse_grid_shape}.")

        squeeze_axis = settings["squeeze_axis"].GetInt()
        if squeeze_axis == -1:
            self.squeeze_axis = None
        elif squeeze_axis in (0, 1, 2):
            self.squeeze_axis = squeeze_axis
        else:
            raise ValueError(f"\"squeeze_axis\" must be -1 (off), 0, 1 or 2, got {squeeze_axis}.")

        self.model_interface = settings["model_interface"].GetString()
        if self.model_interface not in ("grid", "modafno"):
            raise ValueError(
                f"Unsupported model interface \"{self.model_interface}\". "
                "Use \"grid\" or \"modafno\".")

        box = settings["bounding_box"].GetVector()
        if len(box) == 0:
            # Default: the fine part's box, so both grids share the domain.
            self.bounding_box = grid_bridge.ComputeBoundingBox(self.fine_model_part)
        elif len(box) == 6:
            self.bounding_box = (numpy.array(box[:3]), numpy.array(box[3:]))
        else:
            raise ValueError("\"bounding_box\" must be empty or [x0,y0,z0,x1,y1,z1].")

        self.execution_point = settings["execution_point"].GetString()
        if self.execution_point not in _EXECUTION_POINTS:
            raise ValueError(
                f"Unsupported execution point \"{self.execution_point}\". "
                f"Supported: {', '.join(_EXECUTION_POINTS)}.")
        self.output_interval = settings["output_interval"].GetInt()
        if self.output_interval < 1:
            raise ValueError(f"\"output_interval\" must be >= 1 [ output_interval = {self.output_interval} ].")

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

    def _ExecuteIfDue(self) -> None:
        if self.coarse_model_part.ProcessInfo[Kratos.STEP] % self.output_interval == 0:
            self.RunSuperResolution()

    def RunSuperResolution(self) -> None:
        torch = torch_bridge._TryImportTorch()
        if self._model is None:
            self._model, self._device = model_registry.LoadModelWithCardCheck(
                self.model_settings, self.input_specs, self.output_specs, type(self).__name__)
            self._normalization = model_registry.LoadOutputNormalization(self.model_settings)
            self._input_normalization = model_registry.LoadInputNormalization(self.model_settings)

        with NvtxRange("PhysicsNeMo::SampleFieldsOnGrid"):
            coarse_grid, _ = grid_bridge.SampleFieldsOnGrid(
                self.coarse_model_part, self.input_specs, self.coarse_grid_shape, self.bounding_box)
        if self.squeeze_axis is not None:
            coarse_grid = coarse_grid.mean(axis=1 + self.squeeze_axis)  # (C, A, B)
        # the card's "input_normalization", per channel along axis 0
        coarse_grid = model_registry.ApplyInputNormalization(
            coarse_grid, self._input_normalization, channel_axis=0)

        if self._ood_guard.enabled:  # grid (C, *spatial) -> (prod(spatial), C) features
            self._ood_guard.Check(
                torch.from_numpy(coarse_grid.reshape(coarse_grid.shape[0], -1).T),
                type(self).__name__)

        with torch.no_grad(), NvtxRange("PhysicsNeMo::Forward"):
            batch = torch.from_numpy(coarse_grid[None]).to(self._device)
            # Match the model's parameter dtype (e.g. float32 conv weights);
            # parameterless models (plain resamplers) keep float64 precision.
            parameter = next(self._model.parameters(), None)
            if parameter is not None:
                batch = batch.to(parameter.dtype)
            if self.model_interface == "modafno":
                mod = torch.full((1, 1), self.coarse_model_part.ProcessInfo[Kratos.TIME],
                                 dtype=batch.dtype, device=batch.device)
                prediction = self._model(batch, mod).cpu().numpy()[0]
            else:
                prediction = self._model(batch).cpu().numpy()[0]

        expected_ndim = 3 if self.squeeze_axis is not None else 4
        if prediction.ndim != expected_ndim:
            layout = "(C, A, B)" if self.squeeze_axis is not None else "(C, D, H, W)"
            raise ValueError(
                f"The model must return a {layout} grid per sample; got shape {prediction.shape}.")
        coarse_sizes = (tuple(n for axis, n in enumerate(self.coarse_grid_shape)
                              if axis != self.squeeze_axis)
                        if self.squeeze_axis is not None else self.coarse_grid_shape)
        for axis, coarse_size in enumerate(coarse_sizes):
            if prediction.shape[1 + axis] % coarse_size != 0:
                Kratos.Logger.PrintWarning(
                    "SuperResolutionProcess",
                    f"Output axis {axis} size {prediction.shape[1 + axis]} is not an integer "
                    f"multiple of the coarse size {coarse_size}.")

        if self.squeeze_axis is not None:
            # duplicate the prediction across the collapsed thin axis
            thin_size = self.coarse_grid_shape[self.squeeze_axis]
            prediction = numpy.repeat(
                numpy.expand_dims(prediction, 1 + self.squeeze_axis), thin_size,
                axis=1 + self.squeeze_axis)

        # the card's "output_normalization" makes a normalized prediction
        # physical - per channel along axis 0, the grid layout
        prediction = model_registry.ApplyOutputNormalization(
            prediction.astype(numpy.float64), self._normalization, channel_axis=0)

        with NvtxRange("PhysicsNeMo::ScatterGridToNodes"):
            grid_bridge.ScatterGridToNodes(
                prediction, self.bounding_box, self.fine_model_part, self.output_specs)
