"""Process writing a time-sampled OpenUSD stage from a running solve.

The digital-twin export: per configured interval the model part's surface
(its tessellation's outward boundary for volume meshes, the triangles
themselves for surface parts) plus the requested NODAL fields are written as
one more time sample of a UsdGeomMesh - a scrubbable, deforming, field-
carrying asset any USD viewer (usdview, Omniverse, Blender, ...) opens
directly. Deployed-surrogate predictions written into ordinary Kratos
variables (InferenceProcess and friends) are exported like any other field,
which is the point: the twin shows what the surrogate says, with its
uncertainty fields alongside if exported too.

"kind": "points" writes a UsdGeomPoints cloud instead (meshless parts -
particle surrogates), gathering node positions and fields directly with no
tessellation, so it also needs neither torch nor physicsnemo.

MPI-aware exactly as CuratorExportProcess: topology and fields are gathered
onto a rank-0 shadow part and rank 0 writes.

The stage is saved in ExecuteFinalize (one file, many time samples - unlike
the per-step files of the curator/vtu exporters). usd-core is imported
lazily at first execution (inside deployment.usd_export).
"""

from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
from KratosMultiphysics.PhysicsNeMoApplication.deployment import usd_export

_SUPPORTED_KINDS = ("mesh", "points")
_TIME_SOURCES = ("step", "time")


def Factory(settings: Kratos.Parameters, model: Kratos.Model) -> "UsdExportProcess":
    if not isinstance(settings, Kratos.Parameters):
        raise TypeError("Expected input shall be a Parameters object, encapsulating a json string")
    return UsdExportProcess(model, settings["Parameters"])


class UsdExportProcess(Kratos.Process):
    """Appends a USD time sample of the model part every output_interval steps."""

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
            "kind"                     : "mesh",
            "source_container"         : "Elements",
            "tessellation_mode"        : "smallest_id_diagonal",
            "higher_order_mode"        : "reduce",
            "curved_refinement_levels" : 2,
            "output_file"              : "physics_nemo_twin.usda",
            "prim_path"                : "",
            "up_axis"                  : "Z",
            "meters_per_unit"          : 1.0,
            "time_source"              : "step",
            "time_codes_per_second"    : 1.0,
            "output_interval"          : 1
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
        for variable_name, data_location in self.field_specs:
            if not data_location.startswith("node_"):
                raise ValueError(
                    f"USD export writes \"vertex\"-interpolated primvars, so only nodal "
                    f"data locations are supported [ field \"{variable_name}\" has "
                    f"data_location \"{data_location}\" ]. Extrapolate element/Gauss "
                    "data to nodes first, or use the vtu/curator exporters.")

        self.kind = settings["kind"].GetString()
        if self.kind not in _SUPPORTED_KINDS:
            raise ValueError(
                f"Unsupported \"kind\" \"{self.kind}\". Supported: {', '.join(_SUPPORTED_KINDS)}.")
        self.source_container = settings["source_container"].GetString()
        self.tessellation_mode = settings["tessellation_mode"].GetString()
        self.higher_order_mode = settings["higher_order_mode"].GetString()
        self.curved_refinement_levels = settings["curved_refinement_levels"].GetInt()
        if self.curved_refinement_levels < 1:
            raise ValueError(
                f"\"curved_refinement_levels\" must be >= 1 "
                f"[ curved_refinement_levels = {self.curved_refinement_levels} ].")

        self.output_file = Path(settings["output_file"].GetString())
        prim_path = settings["prim_path"].GetString()
        self.prim_path = prim_path or f"/Kratos/{self.model_part.Name}"
        if not self.prim_path.startswith("/"):
            raise ValueError(f"\"prim_path\" must be absolute, got \"{self.prim_path}\".")
        self.up_axis = settings["up_axis"].GetString()
        self.meters_per_unit = settings["meters_per_unit"].GetDouble()
        self.time_source = settings["time_source"].GetString()
        if self.time_source not in _TIME_SOURCES:
            raise ValueError(
                f"Unsupported \"time_source\" \"{self.time_source}\". "
                f"Supported: {', '.join(_TIME_SOURCES)}.")
        self.time_codes_per_second = settings["time_codes_per_second"].GetDouble()
        self.output_interval = settings["output_interval"].GetInt()
        if self.output_interval < 1:
            raise ValueError(f"\"output_interval\" must be >= 1 [ output_interval = {self.output_interval} ].")

        self._stage = None
        self._provenance_cache = None

    def _GetStage(self):
        # Built on first use, so constructing the process never requires
        # usd-core (the lazy-import contract).
        if self._stage is None:
            self._stage = usd_export.CreateUsdStage(
                self.output_file, up_axis=self.up_axis,
                meters_per_unit=self.meters_per_unit,
                time_codes_per_second=self.time_codes_per_second)
        return self._stage

    def ExecuteInitialize(self) -> None:
        if self.output_file.parent != Path("."):
            self.output_file.parent.mkdir(parents=True, exist_ok=True)

    def ExecuteFinalizeSolutionStep(self) -> None:
        step = self.model_part.ProcessInfo[Kratos.STEP]
        if step % self.output_interval != 0:
            return

        if self.model_part.IsDistributed():
            # Same contract as the other exporters: all ranks enter the
            # collectives, only rank 0 comes back with a part to write.
            from KratosMultiphysics.PhysicsNeMoApplication.distributed import distributed_utils
            gathered = distributed_utils.GatherModelPartToRank0(
                self.model_part, self.field_specs, self.source_container)
            if gathered.model_part is None:
                return
            export_part, export_field_specs = gathered.model_part, gathered.field_specs
        else:
            export_part, export_field_specs = self.model_part, self.field_specs

        time_code = (float(step) if self.time_source == "step"
                     else float(self.model_part.ProcessInfo[Kratos.TIME]))
        if self.kind == "points":
            self._WritePointsSample(export_part, export_field_specs, time_code)
        else:
            self._WriteMeshSample(export_part, export_field_specs, time_code)

    def ExecuteFinalize(self) -> None:
        self.Save()

    def Save(self) -> None:
        """Persists the stage; a no-op before the first written sample."""
        if self._stage is not None:
            usd_export.SaveStage(self._stage)

    def _WritePointsSample(self, export_part, field_specs, time_code) -> None:
        from KratosMultiphysics.PhysicsNeMoApplication.training.streaming_dataset import GatherSampleArrays

        points = numpy.array([[node.X, node.Y, node.Z] for node in export_part.Nodes])
        arrays = GatherSampleArrays(export_part, field_specs)
        point_fields = {
            variable_name: arrays[f"{variable_name}__{data_location}"]
            for variable_name, data_location in field_specs
        }
        usd_export.WritePointsTimeSample(
            self._GetStage(), self.prim_path, points, point_fields, time_code)

    def _WriteMeshSample(self, export_part, field_specs, time_code) -> None:
        from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import domain_mesh_builder
        from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import spatial

        if self._provenance_cache is None:
            self._provenance_cache = domain_mesh_builder.ProvenanceCache(
                self.source_container, self.tessellation_mode,
                self.higher_order_mode, self.curved_refinement_levels)
        resolved_specs = [
            (Kratos.KratosGlobals.GetVariable(variable_name), data_location)
            for variable_name, data_location in field_specs
        ]
        cached = (self._provenance_cache.Get(export_part)
                  if not self.model_part.IsDistributed() else None)
        mesh, _ = domain_mesh_builder.BuildMesh(
            export_part, resolved_specs, self.source_container,
            self.tessellation_mode, self.higher_order_mode, self.curved_refinement_levels,
            provenance=cached)

        points = mesh.points.detach().cpu().numpy()
        cell_width = int(mesh.cells.shape[1])
        if cell_width == 4:
            # volume mesh: its outward-oriented boundary is the visible
            # surface; ALL points are kept so field arrays stay aligned
            triangles = spatial.BoundarySurface(mesh).cells.detach().cpu().numpy()
        elif cell_width == 3:
            triangles = mesh.cells.detach().cpu().numpy()
        else:
            raise ValueError(
                f"Unexpected tessellation cell width {cell_width}; the mesh bridge "
                "produces triangles or tetrahedra.")
        point_fields = {
            variable.Name(): mesh.point_data[variable.Name()].detach().cpu().numpy()
            for variable, _ in resolved_specs
        }
        usd_export.WriteMeshTimeSample(
            self._GetStage(), self.prim_path, points, triangles, point_fields, time_code)
