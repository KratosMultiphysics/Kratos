"""Process deploying same-resolution grid-to-grid models (FNO, UNet, ...).

The pattern of PhysicsNeMo's Darcy-FNO and datacenter-thermal-UNet examples:
sample the input fields onto a regular (C, D, H, W) grid over the model part,
run a model that returns a grid of the SAME spatial size, and scatter the
output fields back onto the same model part's nodes.

Implemented as the single-model-part special case of SuperResolutionProcess
(same source and target part, scaling factor 1) — all sampling, dtype and
scatter machinery is shared.

torch/physicsnemo are imported lazily on first execution.
"""

import KratosMultiphysics as Kratos
from KratosMultiphysics.PhysicsNeMoApplication.processes.inference.superresolution_process import SuperResolutionProcess


def Factory(settings: Kratos.Parameters, model: Kratos.Model) -> "GridInferenceProcess":
    if not isinstance(settings, Kratos.Parameters):
        raise TypeError("Expected input shall be a Parameters object, encapsulating a json string")
    return GridInferenceProcess(model, settings["Parameters"])


class GridInferenceProcess(SuperResolutionProcess):
    """Grid-in/grid-out inference on one model part.

    Settings are those of SuperResolutionProcess with "model_part_name" and
    "grid_shape" replacing the coarse/fine pair:

        {
            "model_part_name" : "...",
            "model_settings"  : { "checkpoint_file" : "fno.mdlus", "checkpoint_type" : "physicsnemo" },
            "input_fields"    : [ ... ],
            "output_fields"   : [ ... ],
            "grid_shape"      : [16, 16, 16],
            "bounding_box"    : [],
            "execution_point" : "finalize_solution_step",
            "output_interval" : 1
        }
    """

    def __init__(self, model: Kratos.Model, settings: Kratos.Parameters) -> None:
        settings = settings.Clone()
        for forbidden in ("coarse_model_part_name", "fine_model_part_name", "coarse_grid_shape"):
            if settings.Has(forbidden):
                raise ValueError(
                    f"\"{forbidden}\" is not a valid setting of GridInferenceProcess; use "
                    "\"model_part_name\"/\"grid_shape\" (or SuperResolutionProcess for the "
                    "two-model-part case).")
        if not settings.Has("model_part_name"):
            raise ValueError("\"model_part_name\" must be specified.")

        model_part_name = settings["model_part_name"].GetString()
        settings.RemoveValue("model_part_name")
        settings.AddString("coarse_model_part_name", model_part_name)
        settings.AddString("fine_model_part_name", model_part_name)
        if settings.Has("grid_shape"):
            grid_shape = settings["grid_shape"].GetVector()
            settings.RemoveValue("grid_shape")
            settings.AddEmptyArray("coarse_grid_shape")
            for entry in grid_shape:
                settings["coarse_grid_shape"].Append(int(entry))
        super().__init__(model, settings)
