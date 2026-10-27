"""Process exporting a series of tessellated meshes through physicsnemo-curator.

Per configured interval, converts the model part (with the requested fields
attached) into a physicsnemo.mesh.Mesh and writes it through one of
curator's mesh sinks: a Zarr store ("<output_path>/<prefix>_<step>.zarr") or
a VTK unstructured grid ("<output_path>/<prefix>_<step>.vtu"). The sinks
name their output from the index they are given, which is the Kratos step,
so a solve produces one store/file per exported step.

This is the in-loop counterpart of curator_bridge.CreateKratosMeshSource:
use the source when curating an existing collection of model parts through
a full Source -> Filter -> Sink pipeline, and this process when a running
solve should emit AI-ready data directly.

MPI-aware in the same way as MeshExportProcess: on a distributed model part
the topology and fields are gathered onto a serial shadow part
(distributed_utils.GatherModelPartToRank0) and rank 0 writes the file with
the exact serial layout.

physicsnemo-curator is imported lazily at first execution (inside the
curator_bridge sink factories).
"""

from pathlib import Path

import KratosMultiphysics as Kratos
from KratosMultiphysics.PhysicsNeMoApplication.bridges import curator_bridge
from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import domain_mesh_builder

_SUPPORTED_SINKS = ("zarr", "vtu")


def Factory(settings: Kratos.Parameters, model: Kratos.Model) -> "CuratorExportProcess":
    if not isinstance(settings, Kratos.Parameters):
        raise TypeError("Expected input shall be a Parameters object, encapsulating a json string")
    return CuratorExportProcess(model, settings["Parameters"])


class CuratorExportProcess(Kratos.Process):
    """Exports the tessellated model part to a curator sink every output_interval steps."""

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
            "sink"              : "zarr",
            "zarr_backend"      : "physicsnemo",
            "output_path"       : "physics_nemo_curated",
            "file_prefix"       : "mesh",
            "output_interval"   : 1,
            "compression_level" : 3,
            "chunk_size_mb"     : 1.0
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
        self.sink_type = settings["sink"].GetString()
        self.zarr_backend = settings["zarr_backend"].GetString()
        if self.zarr_backend not in ("physicsnemo", "curator"):
            raise ValueError(
                f"Unsupported zarr_backend \"{self.zarr_backend}\". Use "
                "\"physicsnemo\" (mesh.io.to_zarr, no curator install needed) or "
                "\"curator\" (physicsnemo-curator's ETL sink).")
        if self.sink_type not in _SUPPORTED_SINKS:
            raise ValueError(
                f"Unsupported \"sink\" \"{self.sink_type}\". "
                f"Supported: {', '.join(_SUPPORTED_SINKS)}.")
        self.output_path = Path(settings["output_path"].GetString())
        self.file_prefix = settings["file_prefix"].GetString()
        self.output_interval = settings["output_interval"].GetInt()
        if self.output_interval < 1:
            raise ValueError(f"\"output_interval\" must be >= 1 [ output_interval = {self.output_interval} ].")
        self.compression_level = settings["compression_level"].GetInt()
        self.chunk_size_mb = settings["chunk_size_mb"].GetDouble()
        self._provenance_cache = domain_mesh_builder.ProvenanceCache(
            self.source_container, self.tessellation_mode,
            self.higher_order_mode, self.curved_refinement_levels)
        self._sink = None

    def _GetSink(self):
        # Built on first use, so constructing the process never requires
        # physicsnemo-curator (the lazy-import contract).
        if self._sink is None:
            naming_template = self.file_prefix + "_{index:04d}"
            if self.sink_type == "zarr":
                self._sink = curator_bridge.CreateZarrSink(
                    self.output_path, compression_level=self.compression_level,
                    chunk_size_mb=self.chunk_size_mb, naming_template=naming_template)
            else:
                self._sink = curator_bridge.CreateVtuSink(
                    self.output_path, naming_template=naming_template)
        return self._sink

    def _WriteZarrDirectly(self, mesh, step: int) -> str:
        """physicsnemo's own Zarr writer, no curator involved.

        physicsnemo-curator is a git-only package that builds an 890 MB Rust
        toolchain at install time, and its sinks are Mesh-typed so a
        DomainMesh fails in them. Since 2.2 the mesh module writes Zarr
        itself, which is what this path uses - and the default, since
        nothing here needs the ETL machinery around it.
        """
        from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import (
            domain_mesh_builder)

        store = self.output_path / (self.file_prefix + f"_{step:04d}.zarr")
        return domain_mesh_builder.SaveMeshZarr(
            mesh, store, zstd_level=self.compression_level)

    def ExecuteInitialize(self) -> None:
        self.output_path.mkdir(parents=True, exist_ok=True)

    def ExecuteFinalizeSolutionStep(self) -> None:
        step = self.model_part.ProcessInfo[Kratos.STEP]
        if step % self.output_interval != 0:
            return

        if self.model_part.IsDistributed():
            # Same contract as MeshExportProcess: all ranks enter the
            # collectives, only rank 0 comes back with a part to write.
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
        cached = (self._provenance_cache.Get(export_part)
                  if not self.model_part.IsDistributed() else None)
        mesh, _ = domain_mesh_builder.BuildMesh(
            export_part, field_specs, self.source_container,
            self.tessellation_mode, self.higher_order_mode, self.curved_refinement_levels,
            provenance=cached)
        # The sink names its output from the index, so passing the step
        # keeps one store/file per exported step instead of overwriting.
        if self.sink_type == "zarr" and self.zarr_backend == "physicsnemo":
            self._WriteZarrDirectly(mesh, step)
        else:
            curator_bridge.WriteMeshToCuratorSink(self._GetSink(), mesh, step)
