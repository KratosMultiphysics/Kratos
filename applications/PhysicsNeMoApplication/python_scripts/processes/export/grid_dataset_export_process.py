"""Process exporting a time series of regular voxel grids for training.

The grid counterpart of DatasetExportProcess: per configured interval, the
requested nodal fields are sampled onto a regular (C, D, H, W) grid through
grid_bridge.SampleFieldsOnGrid and written to
"<output_path>/<file_prefix>_<step>.npz" with keys "grid" (float32),
"TIME", "STEP" and "bounding_box". The bounding box is resolved once at the
first export and reused for every later step, so all grids of a series live
on the same lattice - the invariant CreateGridSequenceDataset and the
grid-sequence models rely on.

For planar 2D cases use the thin-axis idiom: a grid_shape like [16, 16, 2]
with a bounding_box slightly padded across the mesh plane.

Pure Kratos + numpy: this module never imports torch or physicsnemo.
"""

from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
from KratosMultiphysics.PhysicsNeMoApplication.bridges import grid_bridge
def Factory(settings: Kratos.Parameters, model: Kratos.Model) -> "GridDatasetExportProcess":
    if not isinstance(settings, Kratos.Parameters):
        raise TypeError("Expected input shall be a Parameters object, encapsulating a json string")
    return GridDatasetExportProcess(model, settings["Parameters"])


class GridDatasetExportProcess(Kratos.Process):
    """Exports the sampled field grid every output_interval steps."""

    def __init__(self, model: Kratos.Model, settings: Kratos.Parameters) -> None:
        super().__init__()

        default_settings = Kratos.Parameters("""{
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "list_of_fields"  : [
                {
                    "variable_name" : "PLEASE_SPECIFY_VARIABLE_NAME",
                    "data_location" : "node_historical"
                }
            ],
            "grid_shape"      : [8, 8, 8],
            "bounding_box"    : [],
            "output_path"     : "physics_nemo_grid_dataset",
            "file_prefix"     : "grid",
            "output_interval" : 1
        }""")
        settings.ValidateAndAssignDefaults(default_settings)
        for i in range(settings["list_of_fields"].size()):
            settings["list_of_fields"][i].ValidateAndAssignDefaults(default_settings["list_of_fields"][0])

        self.model_part = model[settings["model_part_name"].GetString()]
        self.field_specs = [
            (settings["list_of_fields"][i]["variable_name"].GetString(),
             settings["list_of_fields"][i]["data_location"].GetString())
            for i in range(settings["list_of_fields"].size())
        ]
        self.grid_shape = tuple(int(n) for n in settings["grid_shape"].GetVector())
        if len(self.grid_shape) != 3:
            raise ValueError(f"\"grid_shape\" must have three entries, got {self.grid_shape}.")

        box = settings["bounding_box"].GetVector()
        if len(box) == 0:
            self.bounding_box = None  # resolved at the first export, then frozen
        elif len(box) == 6:
            self.bounding_box = (numpy.array(box[:3]), numpy.array(box[3:]))
        else:
            raise ValueError("\"bounding_box\" must be empty or [x0,y0,z0,x1,y1,z1].")

        self.output_path = Path(settings["output_path"].GetString())
        self.file_prefix = settings["file_prefix"].GetString()
        self.output_interval = settings["output_interval"].GetInt()
        if self.output_interval < 1:
            raise ValueError(f"\"output_interval\" must be >= 1 [ output_interval = {self.output_interval} ].")

    def ExecuteInitialize(self) -> None:
        self.output_path.mkdir(parents=True, exist_ok=True)

    def ExecuteFinalizeSolutionStep(self) -> None:
        step = self.model_part.ProcessInfo[Kratos.STEP]
        if step % self.output_interval != 0:
            return

        if self.bounding_box is None:
            self.bounding_box = grid_bridge.ComputeBoundingBox(self.model_part)

        grid, _ = grid_bridge.SampleFieldsOnGrid(
            self.model_part, self.field_specs, self.grid_shape, self.bounding_box)
        numpy.savez(
            self.output_path / f"{self.file_prefix}_{step}.npz",
            grid=grid.astype(numpy.float32),
            TIME=self.model_part.ProcessInfo[Kratos.TIME],
            STEP=step,
            bounding_box=numpy.concatenate([numpy.asarray(b, dtype=float) for b in self.bounding_box]))
