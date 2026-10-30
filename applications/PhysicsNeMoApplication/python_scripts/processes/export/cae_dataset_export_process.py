"""Process exporting cases in physicsnemo's CAE datapipe layout.

Writes one ``.npz`` per case in the per-sample dictionary format
``physicsnemo.datapipes.cae.CAEDataset`` reads, carrying the superset of the
keys ``DoMINODataPipe`` and ``TransolverDataPipe`` consume (each pipe
selects its own subset through ``keys_to_read``, extra keys are ignored):

- STL geometry from the tessellated surface: ``stl_coordinates`` (P, 3),
  ``stl_faces`` (3T,) int32 (FLATTENED - DoMINO expects it flat),
  ``stl_centers`` (T, 3), ``stl_areas`` (T,), plus the surface aliases
  ``surface_mesh_centers`` / ``surface_normals`` / ``surface_areas``.
- ``surface_fields`` (T, F_s) / ``volume_fields`` (N, F_v) float32 (keys
  omitted when the corresponding field list is empty - DoMINO then runs in
  inference mode).
- ``volume_mesh_centers`` (N, 3): the volume part's nodes.
- Global parameters: every key of ``global_params`` is written BOTH as an
  individual shape-(1,) float32 array (Transolver style, e.g.
  ``stream_velocity``/``air_density``) AND stacked into
  ``global_params_values`` / ``global_params_reference`` (k, 1) float32
  (DoMINO style). The stacking order is ``global_params_order`` when given
  (recommended for DoMINO, which is order-sensitive - e.g. velocity first),
  alphabetical otherwise (Kratos Parameters do not preserve insertion
  order).
- ``TIME`` / ``STEP`` as shape-(1,) arrays (0-d arrays break the reader).

The surface triangulation reuses the mesh bridge (BuildProvenance), so the
watertight smallest-id tessellation, higher-order handling and provenance
all apply. MPI-aware: fields and topology are gathered and rank 0 writes
the file with the exact serial layout.

Pure Kratos + numpy: this module never imports torch or physicsnemo. The
matching datapipe factories live in torch_dataset (CreateCaeDataset,
CreateDoMINODataPipe, CreateTransolverDataPipe).
"""

from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import domain_mesh_builder
from KratosMultiphysics.PhysicsNeMoApplication.utilities.tensor_adaptor_dataset_utils import GetTensorAdaptor

_NODAL_LOCATIONS = ("node_historical", "node_non_historical")


def Factory(settings: Kratos.Parameters, model: Kratos.Model) -> "CaeDatasetExportProcess":
    if not isinstance(settings, Kratos.Parameters):
        raise TypeError("Expected input shall be a Parameters object, encapsulating a json string")
    return CaeDatasetExportProcess(model, settings["Parameters"])


def _ReadFreeFormDoubleDict(settings: Kratos.Parameters, key: str) -> dict:
    """Extracts a free-form {name: double} sub-block (insertion-ordered)
    before the fixed-schema validation runs."""
    values = {}
    if settings.Has(key):
        block = settings[key]
        for name in block.keys():
            values[name] = block[name].GetDouble()
        settings.RemoveValue(key)
    return values


class CaeDatasetExportProcess(Kratos.Process):
    """Exports CAE datapipe cases every output_interval steps."""

    def __init__(self, model: Kratos.Model, settings: Kratos.Parameters) -> None:
        super().__init__()

        # free-form dicts leave the settings before the schema validation
        self.global_params = _ReadFreeFormDoubleDict(settings, "global_params")
        self.global_params_reference = _ReadFreeFormDoubleDict(settings, "global_params_reference")
        for name in self.global_params_reference:
            if name not in self.global_params:
                raise ValueError(
                    f"global_params_reference key \"{name}\" has no matching global_params entry.")

        default_settings = Kratos.Parameters("""{
            "model_part_name"          : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "surface_model_part_name"  : "PLEASE_SPECIFY_SURFACE_MODEL_PART_NAME",
            "surface_source_container" : "Conditions",
            "surface_fields"           : [],
            "volume_fields"            : [],
            "global_params_order"      : [],
            "tessellation_mode"        : "smallest_id_diagonal",
            "higher_order_mode"        : "reduce",
            "curved_refinement_levels" : 2,
            "output_path"              : "physics_nemo_cae_dataset",
            "file_prefix"              : "case",
            "case_id"                  : -1,
            "output_interval"          : 1
        }""")
        settings.ValidateAndAssignDefaults(default_settings)
        field_defaults = Kratos.Parameters("""{
            "variable_name" : "PLEASE_SPECIFY_VARIABLE_NAME",
            "data_location" : "node_historical"
        }""")
        for key in ("surface_fields", "volume_fields"):
            for i in range(settings[key].size()):
                settings[key][i].ValidateAndAssignDefaults(field_defaults)

        self.model_part = model[settings["model_part_name"].GetString()]
        self.surface_model_part = model[settings["surface_model_part_name"].GetString()]
        self.surface_source_container = settings["surface_source_container"].GetString()
        if self.surface_source_container not in ("Elements", "Conditions"):
            raise ValueError(
                f"Unsupported surface source container \"{self.surface_source_container}\". "
                "Use \"Elements\" or \"Conditions\".")

        def read_specs(key):
            return [(settings[key][i]["variable_name"].GetString(),
                     settings[key][i]["data_location"].GetString())
                    for i in range(settings[key].size())]

        self.surface_field_specs = read_specs("surface_fields")
        self.volume_field_specs = read_specs("volume_fields")
        for variable_name, data_location in self.volume_field_specs:
            if data_location not in _NODAL_LOCATIONS:
                raise ValueError(
                    f"Volume field \"{variable_name}\" has location \"{data_location}\"; volume "
                    f"centers are the nodes, so only {_NODAL_LOCATIONS} are supported.")

        self.global_params_order = settings["global_params_order"].GetStringArray()
        if self.global_params_order:
            if sorted(self.global_params_order) != sorted(self.global_params):
                raise ValueError(
                    f"\"global_params_order\" {self.global_params_order} must name exactly the "
                    f"global_params keys {sorted(self.global_params)}.")
        else:
            self.global_params_order = sorted(self.global_params)

        self.tessellation_mode = settings["tessellation_mode"].GetString()
        self.higher_order_mode = settings["higher_order_mode"].GetString()
        self.curved_refinement_levels = settings["curved_refinement_levels"].GetInt()
        self.output_path = Path(settings["output_path"].GetString())
        self.file_prefix = settings["file_prefix"].GetString()
        self.case_id = settings["case_id"].GetInt()
        self.output_interval = settings["output_interval"].GetInt()
        self._provenance_cache = domain_mesh_builder.ProvenanceCache(
            self.surface_source_container, self.tessellation_mode,
            self.higher_order_mode, self.curved_refinement_levels)
        self._entity_rows = None
        self._entity_rows_for = None
        if self.output_interval < 1:
            raise ValueError(f"\"output_interval\" must be >= 1 [ output_interval = {self.output_interval} ].")

    def ExecuteInitialize(self) -> None:
        self.output_path.mkdir(parents=True, exist_ok=True)

    def ExecuteFinalizeSolutionStep(self) -> None:
        step = self.model_part.ProcessInfo[Kratos.STEP]
        if step % self.output_interval != 0:
            return
        self.Export()

    # --- serial building blocks ---------------------------------------------

    @staticmethod
    def _SurfaceArrays(provenance):
        if provenance.simplex_cells.shape[1] != 3:
            raise RuntimeError(
                "The surface part must tessellate into triangles (surface geometries); got "
                f"simplices with {provenance.simplex_cells.shape[1]} nodes - check "
                "surface_model_part_name/surface_source_container.")
        triangles = provenance.simplex_points[provenance.simplex_cells]  # (T, 3, 3)
        cross = numpy.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0])
        areas = 0.5 * numpy.linalg.norm(cross, axis=-1)
        normals = cross / (2.0 * areas)[:, None]
        return {
            "stl_coordinates": provenance.simplex_points.astype(numpy.float32),
            "stl_faces": provenance.simplex_cells.astype(numpy.int32).ravel(),  # FLATTENED
            "stl_centers": triangles.mean(axis=1).astype(numpy.float32),
            "stl_areas": areas.astype(numpy.float32),
            "surface_normals": normals.astype(numpy.float32),
        }

    @staticmethod
    def _SurfaceFieldBlock(surface_part, source_container, provenance, field_specs,
                           entity_rows=None):
        """(T, sum widths) float32: per-triangle values of the surface fields."""
        entity_location = "condition" if source_container == "Conditions" else "element"
        node_ids = [node.Id for node in surface_part.Nodes]
        blocks = []
        for variable_name, data_location in field_specs:
            variable = Kratos.KratosGlobals.GetVariable(variable_name)
            tensor_adaptor = GetTensorAdaptor(surface_part, data_location, variable)
            data = numpy.array(tensor_adaptor.data, dtype=numpy.float64)
            if data_location in _NODAL_LOCATIONS:
                gathered = provenance.GatherNodalField(node_ids, data.reshape(len(node_ids), -1))
                per_triangle = gathered[provenance.simplex_cells].mean(axis=1)  # vertex mean
            else:
                if data_location.endswith("_gauss_point"):
                    data = data.mean(axis=1)  # collapse gauss points first
                data = data.reshape(data.shape[0], -1)
                if entity_rows is None:
                    container = (surface_part.Conditions if source_container == "Conditions"
                                 else surface_part.Elements)
                    row_of_entity = {entity.Id: row for row, entity in enumerate(container)}
                    entity_rows = numpy.fromiter(
                        (row_of_entity[int(entity_id)] for entity_id in provenance.cell_provenance[:, 0]),
                        dtype=numpy.int64, count=provenance.number_of_cells)
                per_triangle = data[entity_rows]  # entity value replicated on its sub-triangles
            blocks.append(per_triangle)
        return numpy.concatenate(blocks, axis=1).astype(numpy.float32)

    @staticmethod
    def _NodalFieldBlock(model_part, field_specs):
        """(N, sum widths) float32 over the part's nodes."""
        blocks = []
        for variable_name, data_location in field_specs:
            variable = Kratos.KratosGlobals.GetVariable(variable_name)
            tensor_adaptor = GetTensorAdaptor(model_part, data_location, variable)
            data = numpy.array(tensor_adaptor.data, dtype=numpy.float64)
            blocks.append(data.reshape(data.shape[0], -1))
        return numpy.concatenate(blocks, axis=1).astype(numpy.float32)

    @staticmethod
    def _NodeCoordinates(model_part):
        position_ta = Kratos.TensorAdaptors.NodePositionTensorAdaptor(
            model_part.Nodes, Kratos.Configuration.Current)
        position_ta.CollectData()
        return numpy.array(position_ta.data, dtype=numpy.float64).reshape(-1, 3)

    def _GlobalParameterArrays(self):
        arrays = {}
        values = []
        references = []
        for name in self.global_params_order:
            value = self.global_params[name]
            arrays[name] = numpy.array([value], dtype=numpy.float32)  # (1,), never 0-d
            values.append(value)
            references.append(self.global_params_reference.get(name, value))
        if values:
            arrays["global_params_values"] = numpy.array(values, dtype=numpy.float32)[:, None]
            arrays["global_params_reference"] = numpy.array(references, dtype=numpy.float32)[:, None]
        return arrays

    # --- export -------------------------------------------------------------

    def Export(self) -> None:
        step = self.model_part.ProcessInfo[Kratos.STEP]

        if self.model_part.IsDistributed():
            from KratosMultiphysics.PhysicsNeMoApplication.distributed import distributed_utils
            # volume: nodes + nodal fields, plain id-sorted arrays
            volume_mesh = distributed_utils.GatherMeshToRank0(self.model_part, None)
            volume_blocks = []
            for variable_name, data_location in self.volume_field_specs:
                _, gathered = distributed_utils.GatherFieldToRank0(
                    self.model_part, variable_name, data_location)
                if gathered is not None:
                    volume_blocks.append(
                        numpy.asarray(gathered, dtype=numpy.float64).reshape(len(gathered), -1))
            # surface: serial pipeline on the gathered shadow part
            gathered_surface = distributed_utils.GatherModelPartToRank0(
                self.surface_model_part, self.surface_field_specs, self.surface_source_container)
            if gathered_surface.model_part is None:
                return  # non-writing rank; collectives are done
            surface_part = gathered_surface.model_part
            surface_field_specs = gathered_surface.field_specs
            volume_coordinates = volume_mesh.coordinates
            volume_block = (numpy.concatenate(volume_blocks, axis=1).astype(numpy.float32)
                            if volume_blocks else None)
        else:
            surface_part = self.surface_model_part
            surface_field_specs = self.surface_field_specs
            volume_coordinates = self._NodeCoordinates(self.model_part)
            volume_block = (self._NodalFieldBlock(self.model_part, self.volume_field_specs)
                            if self.volume_field_specs else None)

        # tessellation is ~99% of an export and is identical while the mesh
        # is; the cache re-tessellates as soon as a node moves, because the
        # map carries simplex_points. Serial only: the distributed branch
        # above gathers a fresh shadow part every export.
        if self.model_part.IsDistributed():
            provenance = domain_mesh_builder.BuildProvenance(
                surface_part, self.surface_source_container,
                self.tessellation_mode, self.higher_order_mode, self.curved_refinement_levels)
        else:
            provenance = self._provenance_cache.Get(surface_part)
        if provenance is not self._entity_rows_for:
            container = (surface_part.Conditions
                         if self.surface_source_container == "Conditions"
                         else surface_part.Elements)
            row_of_entity = {entity.Id: row for row, entity in enumerate(container)}
            self._entity_rows = numpy.fromiter(
                (row_of_entity[int(entity_id)] for entity_id in provenance.cell_provenance[:, 0]),
                dtype=numpy.int64, count=provenance.number_of_cells)
            self._entity_rows_for = provenance

        arrays = self._SurfaceArrays(provenance)
        arrays["surface_mesh_centers"] = arrays["stl_centers"]
        arrays["surface_areas"] = arrays["stl_areas"]
        if surface_field_specs:
            arrays["surface_fields"] = self._SurfaceFieldBlock(
                surface_part, self.surface_source_container, provenance,
                surface_field_specs, entity_rows=self._entity_rows)
        arrays["volume_mesh_centers"] = volume_coordinates.astype(numpy.float32)
        if volume_block is not None:
            arrays["volume_fields"] = volume_block
        arrays.update(self._GlobalParameterArrays())
        arrays["TIME"] = numpy.array([self.model_part.ProcessInfo[Kratos.TIME]], dtype=numpy.float32)
        arrays["STEP"] = numpy.array([step], dtype=numpy.int64)

        suffix = self.case_id if self.case_id >= 0 else step
        numpy.savez(self.output_path / f"{self.file_prefix}_{suffix}.npz", **arrays)
