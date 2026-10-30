"""Process exporting a series of tessellated meshes for mesh-based training.

Per configured interval, converts the model part (with the requested fields
attached) into a physicsnemo.mesh.Mesh and saves it to
"<output_path>/<file_prefix>_<step>.pmsh" — the ".pmsh" suffix is what
physicsnemo.datapipes.mesh_dataset.MeshReader's default glob looks for, so
the output directory is directly consumable by CreateMeshDataset /
MeshDataset for training mesh-based models.

MPI-aware: on a distributed model part the topology and fields are gathered
onto a serial shadow part (distributed_utils.GatherModelPartToRank0) and
rank 0 writes the file with the exact serial layout, so consumers are
rank-count-agnostic.

physicsnemo is imported lazily at first execution (inside BuildMesh).
"""

from pathlib import Path

import KratosMultiphysics as Kratos
from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import domain_mesh_builder


def Factory(settings: Kratos.Parameters, model: Kratos.Model) -> "MeshExportProcess":
    if not isinstance(settings, Kratos.Parameters):
        raise TypeError("Expected input shall be a Parameters object, encapsulating a json string")
    return MeshExportProcess(model, settings["Parameters"])


class MeshExportProcess(Kratos.Process):
    """Exports the tessellated model part every output_interval steps."""

    def __init__(self, model: Kratos.Model, settings: Kratos.Parameters) -> None:
        super().__init__()

        default_settings = Kratos.Parameters("""{
            "model_part_name"  : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "list_of_fields"   : [
                {
                    "variable_name" : "PLEASE_SPECIFY_VARIABLE_NAME",
                    "data_location" : "node_historical"
                }
            ],
            "source_container"  : "Elements",
            "tessellation_mode"        : "smallest_id_diagonal",
            "higher_order_mode"        : "reduce",
            "curved_refinement_levels" : 2,
            "output_path"       : "physics_nemo_meshes",
            "file_prefix"       : "mesh",
            "output_interval"   : 1,
            "zero_pad_steps"    : 0
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
        self.source_container = settings["source_container"].GetString()
        self.tessellation_mode = settings["tessellation_mode"].GetString()
        self.higher_order_mode = settings["higher_order_mode"].GetString()
        self.curved_refinement_levels = settings["curved_refinement_levels"].GetInt()
        if self.curved_refinement_levels < 1:
            raise ValueError(
                f"\"curved_refinement_levels\" must be >= 1 "
                f"[ curved_refinement_levels = {self.curved_refinement_levels} ].")
        self.output_path = Path(settings["output_path"].GetString())
        self.file_prefix = settings["file_prefix"].GetString()
        self.output_interval = settings["output_interval"].GetInt()
        self.zero_pad_steps = settings["zero_pad_steps"].GetInt()
        self._provenance_cache = domain_mesh_builder.ProvenanceCache(
            self.source_container, self.tessellation_mode,
            self.higher_order_mode, self.curved_refinement_levels)
        if self.zero_pad_steps < 0:
            raise ValueError(f"\"zero_pad_steps\" must be >= 0 [ zero_pad_steps = {self.zero_pad_steps} ].")
        if self.output_interval < 1:
            raise ValueError(f"\"output_interval\" must be >= 1 [ output_interval = {self.output_interval} ].")

    def ExecuteInitialize(self) -> None:
        self.output_path.mkdir(parents=True, exist_ok=True)

    def ExecuteFinalizeSolutionStep(self) -> None:
        step = self.model_part.ProcessInfo[Kratos.STEP]
        if step % self.output_interval != 0:
            return

        if self.model_part.IsDistributed():
            # Gather topology + fields onto a serial "shadow" part; rank 0
            # runs the ordinary serial export on it (id-sorted construction
            # makes the output identical to a serial run's), other ranks
            # only participate in the collectives. The topology is
            # re-gathered every output step, so moving meshes stay correct.
            from KratosMultiphysics.PhysicsNeMoApplication.distributed import distributed_utils
            gathered = distributed_utils.GatherModelPartToRank0(
                self.model_part, self.field_specs, self.source_container)
            if gathered.model_part is None:
                return
            export_part, export_field_specs = gathered.model_part, gathered.field_specs
        else:
            export_part, export_field_specs = self.model_part, self.field_specs

        field_specs = [
            (Kratos.KratosGlobals.GetVariable(variable_name), data_location)
            for variable_name, data_location in export_field_specs
        ]
        # the tessellation is reused while the mesh is unchanged; the field
        # data is always collected fresh, which is the point of the export.
        # Serial only: the distributed branch above hands back a newly
        # gathered shadow part every step.
        cached = (self._provenance_cache.Get(export_part)
                  if not self.model_part.IsDistributed() else None)
        mesh, _ = domain_mesh_builder.BuildMesh(
            export_part, field_specs, self.source_container,
            self.tessellation_mode, self.higher_order_mode, self.curved_refinement_levels,
            provenance=cached)
        # ".pmsh" suffix is required for MeshReader's default "**/*.pmsh" glob.
        # MeshReader sorts its glob lexicographically, so an unpadded series
        # reads back as 1, 10, 11, 2, ... - set "zero_pad_steps" for
        # trajectory exports whose item order must be the time order.
        name = f"{step:0{self.zero_pad_steps}d}" if self.zero_pad_steps else str(step)
        domain_mesh_builder.SaveMesh(mesh, self.output_path / f"{self.file_prefix}_{name}.pmsh")
