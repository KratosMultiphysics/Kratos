"""Builds physicsnemo meshes from Kratos model parts and scatters fields back.

physicsnemo is an optional runtime dependency: everything except
BuildDomainMesh works without it (BuildProvenance and ScatterFieldBack are
pure Kratos + numpy).
"""

import numpy

import KratosMultiphysics as Kratos

from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import tessellation
from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge.provenance import MeshProvenanceMap
from KratosMultiphysics.PhysicsNeMoApplication.utilities.tensor_adaptor_dataset_utils import GetTensorAdaptor


def _TryImportPhysicsNemo():
    try:
        import physicsnemo.mesh
        return physicsnemo
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.mesh_bridge.domain_mesh_builder.BuildDomainMesh requires "
            "physicsnemo, which could not be imported. Install it with e.g. "
            "'pip install nvidia-physicsnemo'.") from e


def _GetContainer(model_part: Kratos.ModelPart, source_container: str):
    if source_container == "Elements":
        return model_part.Elements
    if source_container == "Conditions":
        return model_part.Conditions
    raise ValueError(f"Unsupported source container \"{source_container}\". Use \"Elements\" or \"Conditions\".")


# Homogeneous containers of these types tessellate as the identity on their
# corner nodes, so the whole provenance can be built with C++ adaptors and
# vectorized numpy - no per-entity python tessellation needed. Higher-order
# simplices qualify only in "reduce" mode (corner truncation).
_VECTORIZABLE_LINEAR_SIMPLICES = frozenset((
    Kratos.GeometryData.KratosGeometryType.Kratos_Triangle2D3,
    Kratos.GeometryData.KratosGeometryType.Kratos_Triangle3D3,
    Kratos.GeometryData.KratosGeometryType.Kratos_Tetrahedra3D4,
))
_VECTORIZABLE_HIGHER_ORDER_SIMPLICES = frozenset((
    Kratos.GeometryData.KratosGeometryType.Kratos_Triangle2D6,
    Kratos.GeometryData.KratosGeometryType.Kratos_Triangle3D6,
    Kratos.GeometryData.KratosGeometryType.Kratos_Tetrahedra3D10,
))


def _TryVectorizedSimplexProvenance(model_part, container, coordinates,
                                    higher_order_mode, source_container):
    """Fast path for homogeneous simplex containers; None when not applicable.

    Produces bit-identical results to the general per-entity path (identity
    tessellation, corner truncation for higher-order types in "reduce"
    mode): entity order, sub-cell indices and the sorted point provenance
    all match by construction.
    """
    geometry_codes = numpy.fromiter(
        (int(entity.GetGeometry().GetGeometryType()) for entity in container),
        dtype=numpy.int64, count=len(container))
    unique_codes = numpy.unique(geometry_codes)
    if len(unique_codes) != 1:
        return None
    geometry_type = Kratos.GeometryData.KratosGeometryType(int(unique_codes[0]))
    if geometry_type in _VECTORIZABLE_LINEAR_SIMPLICES:
        pass
    elif geometry_type in _VECTORIZABLE_HIGHER_ORDER_SIMPLICES and higher_order_mode == "reduce":
        pass
    else:
        return None

    corner_count = tessellation._CORNER_COUNT[geometry_type]
    connectivity_ta = Kratos.TensorAdaptors.ConnectivityIdsTensorAdaptor(container)
    connectivity_ta.CollectData()
    corners = numpy.array(connectivity_ta.data, dtype=numpy.int64)[:, :corner_count]

    entity_ids = numpy.fromiter(
        (entity.Id for entity in container), dtype=numpy.int64, count=len(container))
    cell_provenance = numpy.stack(
        [entity_ids, numpy.zeros(len(entity_ids), dtype=numpy.int64)], axis=1)

    part_node_ids = numpy.fromiter(
        (node.Id for node in model_part.Nodes), dtype=numpy.int64,
        count=model_part.NumberOfNodes())
    referenced = numpy.unique(corners)  # sorted, matching the general path
    if len(part_node_ids) == 0 or numpy.any(numpy.diff(part_node_ids) <= 0):
        return None  # unsorted node container: let the general path handle it
    rows = numpy.searchsorted(part_node_ids, referenced)
    if rows[-1] >= len(part_node_ids) or not numpy.array_equal(part_node_ids[rows], referenced):
        return None  # entity references nodes outside the part

    return MeshProvenanceMap(
        simplex_points=numpy.ascontiguousarray(coordinates[rows]),
        simplex_cells=numpy.searchsorted(referenced, corners),
        point_provenance=referenced,
        cell_provenance=cell_provenance,
        source_container=source_container)


def BuildProvenance(model_part: Kratos.ModelPart, source_container: str = "Elements",
                    tessellation_mode: str = "smallest_id_diagonal",
                    higher_order_mode: str = "reduce",
                    curved_refinement_levels: int = 2) -> MeshProvenanceMap:
    """Tessellates a model part into simplices with full field provenance.

    Homogeneous simplex containers (linear triangles/tetrahedra, and their
    quadratic variants in "reduce" mode) take a vectorized fast path built
    on the C++ connectivity adaptor; everything else goes through the
    per-entity tessellation tables. higher_order_mode="curved" samples the
    exact isoparametric geometry with synthetic points (see
    curved_tessellation) at 2^curved_refinement_levels cells per parametric
    axis.

    Args:
        model_part: The model part to tessellate.
        source_container: "Elements" (default) or "Conditions".
        tessellation_mode: "smallest_id_diagonal" (default) or "fan", see
            tessellation.TessellateEntity.
        higher_order_mode: "reduce" (default), "subdivide" or "curved".
        curved_refinement_levels: parametric refinement depth of the curved
            mode (>= 1; ignored by the other modes).

    Returns:
        MeshProvenanceMap for the tessellated mesh.
    """
    container = _GetContainer(model_part, source_container)
    if len(container) == 0:
        raise RuntimeError(
            f"Model part \"{model_part.FullName()}\" has no entities in its {source_container} container.")

    position_ta = Kratos.TensorAdaptors.NodePositionTensorAdaptor(model_part.Nodes, Kratos.Configuration.Current)
    position_ta.CollectData()
    coordinates = numpy.array(position_ta.data)

    if higher_order_mode == "curved":
        if tessellation_mode != "smallest_id_diagonal":
            raise ValueError(
                "higher_order_mode=\"curved\" requires tessellation_mode="
                "\"smallest_id_diagonal\" (its diagonal choices are key-driven).")
        if int(curved_refinement_levels) < 1:
            raise ValueError(
                f"\"curved_refinement_levels\" must be >= 1, got {curved_refinement_levels}.")
    else:
        fast = _TryVectorizedSimplexProvenance(
            model_part, container, coordinates, higher_order_mode, source_container)
        if fast is not None:
            return fast

    node_ids = [node.Id for node in model_part.Nodes]
    node_coordinates = {node_id: coordinates[row] for row, node_id in enumerate(node_ids)}

    if higher_order_mode == "curved":
        from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import curved_tessellation
        result = curved_tessellation.TessellateContainerCurved(
            container, node_coordinates, curved_refinement_levels)
        return MeshProvenanceMap(
            simplex_points=result.point_coordinates,
            simplex_cells=result.simplex_cells,
            point_provenance=result.point_node_ids,
            cell_provenance=result.cell_provenance,
            source_container=source_container,
            synthetic_parent_ids=result.synthetic_parent_ids,
            synthetic_local_coordinates=result.synthetic_local_coordinates,
            synthetic_node_ids=result.synthetic_node_ids,
            synthetic_weights=result.synthetic_weights)

    simplex_node_ids, cell_provenance = tessellation.TessellateContainer(
        container, node_coordinates, tessellation_mode, higher_order_mode)

    # Corner-only tessellation: keep only the nodes actually referenced.
    referenced = sorted({node_id for simplex in simplex_node_ids for node_id in simplex})
    point_provenance = numpy.array(referenced, dtype=numpy.int64)
    node_id_to_point = {node_id: index for index, node_id in enumerate(referenced)}

    simplex_points = numpy.array([node_coordinates[node_id] for node_id in referenced], dtype=numpy.float64)
    simplex_cells = numpy.array(
        [[node_id_to_point[node_id] for node_id in simplex] for simplex in simplex_node_ids],
        dtype=numpy.int64)

    return MeshProvenanceMap(
        simplex_points=simplex_points,
        simplex_cells=simplex_cells,
        point_provenance=point_provenance,
        cell_provenance=cell_provenance,
        source_container=source_container)


def _GatherCellRows(model_part: Kratos.ModelPart, provenance: MeshProvenanceMap) -> numpy.ndarray:
    """Row index (into the source container) of each simplex cell's entity."""
    container = _GetContainer(model_part, provenance.source_container)
    container_ids = numpy.fromiter(
        (entity.Id for entity in container), dtype=numpy.int64, count=len(container))
    source_ids = provenance.cell_provenance[:, 0]
    if len(container_ids) > 0 and numpy.all(numpy.diff(container_ids) > 0):
        positions = numpy.searchsorted(container_ids, source_ids)
        clipped = numpy.minimum(positions, len(container_ids) - 1)
        if not numpy.array_equal(container_ids[clipped], source_ids):
            raise KeyError("Provenance references entities missing from the container.")
        return clipped
    entity_rows = {int(entity.Id): row for row, entity in enumerate(container)}
    return numpy.fromiter(
        (entity_rows[int(entity_id)] for entity_id in source_ids),
        dtype=numpy.int64, count=provenance.number_of_cells)


def _CollectFieldData(model_part: Kratos.ModelPart, field_specs, provenance: MeshProvenanceMap):
    """Gathers (point_data, cell_data) dicts onto the tessellated mesh."""
    node_ids = [node.Id for node in model_part.Nodes]
    point_data = {}
    cell_data = {}
    for variable, data_location in field_specs:
        tensor_adaptor = GetTensorAdaptor(model_part, data_location, variable)
        data = numpy.array(tensor_adaptor.data)
        if data_location in ("node_historical", "node_non_historical"):
            point_data[variable.Name()] = provenance.GatherNodalField(node_ids, data)
        elif data_location in ("element", "condition"):
            cell_data[variable.Name()] = data[_GatherCellRows(model_part, provenance)]
        elif data_location in ("element_gauss_point", "condition_gauss_point"):
            # Gauss points have no geometric counterpart in the simplex mesh:
            # collapse to a per-entity mean (axis 1 = gauss point index).
            cell_data[variable.Name()] = data.mean(axis=1)[_GatherCellRows(model_part, provenance)]
        else:
            raise ValueError(f"Unsupported data location \"{data_location}\".")
    return point_data, cell_data


def _MeshFromProvenance(physicsnemo, provenance: MeshProvenanceMap, point_data, cell_data):
    import torch  # physicsnemo guarantees torch is present
    return physicsnemo.mesh.Mesh(
        points=torch.from_numpy(provenance.simplex_points),
        cells=torch.from_numpy(provenance.simplex_cells),
        point_data={name: torch.as_tensor(value) for name, value in point_data.items()},
        cell_data={name: torch.as_tensor(value) for name, value in cell_data.items()})


class ProvenanceCache:
    """Reuses a provenance map until the mesh it describes changes.

    Tessellation dominates the export processes - `BuildProvenance` measures
    at ~143 ms of a ~145 ms surface export in benchmarks/benchmark_bridges.py
    terms - and the result is identical while the mesh is.

    Invalidation is by **node coordinates**, not just entity count, because a
    `MeshProvenanceMap` carries `simplex_points`: reusing one across a
    deforming mesh would export stale geometry. Gathering and comparing the
    coordinates costs ~0.16 ms against a ~143 ms rebuild, so the check is
    close to free and exact - a 1e-9 move invalidates.

    Instantiate one per process. This is deliberately a class rather than a
    module-level cache: the lifetime of the state belongs to the caller.
    """

    def __init__(self, source_container: str = "Elements",
                 tessellation_mode: str = "smallest_id_diagonal",
                 higher_order_mode: str = "reduce",
                 curved_refinement_levels: int = 2):
        self.source_container = source_container
        self.tessellation_mode = tessellation_mode
        self.higher_order_mode = higher_order_mode
        self.curved_refinement_levels = curved_refinement_levels
        self._provenance = None
        self._coordinates = None
        self._number_of_entities = None

    @staticmethod
    def _Coordinates(model_part: Kratos.ModelPart):
        adaptor = Kratos.TensorAdaptors.NodePositionTensorAdaptor(
            model_part.Nodes, Kratos.Configuration.Current)
        adaptor.CollectData()
        return numpy.array(adaptor.data)

    def Get(self, model_part: Kratos.ModelPart):
        """The provenance map for this model part, tessellating only if needed."""
        coordinates = self._Coordinates(model_part)
        number_of_entities = len(_GetContainer(model_part, self.source_container))
        if (self._provenance is not None
                and self._number_of_entities == number_of_entities
                and numpy.array_equal(self._coordinates, coordinates)):
            return self._provenance

        self._provenance = BuildProvenance(
            model_part, self.source_container, self.tessellation_mode,
            self.higher_order_mode, self.curved_refinement_levels)
        self._coordinates = coordinates
        self._number_of_entities = number_of_entities
        return self._provenance


def BuildMesh(model_part: Kratos.ModelPart, field_specs=(), source_container: str = "Elements",
              tessellation_mode: str = "smallest_id_diagonal", higher_order_mode: str = "reduce",
              curved_refinement_levels: int = 2, provenance=None):
    """Builds a physicsnemo.mesh.Mesh (plus provenance) from a model part.

    Args:
        model_part: The model part to convert.
        field_specs: iterable of (variable, data_location) pairs attached as
            point_data (nodal locations) or cell_data (element locations).
            Gauss-point locations are collapsed to a per-source-entity mean
            and attached as cell_data (constant over each entity's sub-cells).
        source_container: "Elements" (default) or "Conditions".
        tessellation_mode: see BuildProvenance.
        higher_order_mode: see BuildProvenance.
        provenance: Optional prebuilt map (see ProvenanceCache) used instead
            of tessellating again. The field data is still collected here, so
            it stays current; only the tessellation is reused.

    Returns:
        (mesh, provenance_map)
    """
    physicsnemo = _TryImportPhysicsNemo()

    if provenance is None:
        provenance = BuildProvenance(model_part, source_container, tessellation_mode,
                                     higher_order_mode, curved_refinement_levels)
    point_data, cell_data = _CollectFieldData(model_part, field_specs, provenance)
    return _MeshFromProvenance(physicsnemo, provenance, point_data, cell_data), provenance


def BuildDomainMesh(model_part: Kratos.ModelPart,
                    field_specs=(),
                    boundary_sub_model_part_names=(),
                    source_container: str = "Elements",
                    tessellation_mode: str = "smallest_id_diagonal",
                    higher_order_mode: str = "reduce",
                    curved_refinement_levels: int = 2):
    """Builds a physicsnemo.mesh.DomainMesh with named boundary meshes.

    The interior mesh is the tessellation of the model part itself; each
    named sub-model-part becomes a named boundary Mesh tessellated from its
    Conditions container (falling back to its Elements container when it has
    no conditions; skipped with a warning when it has neither). The same
    field_specs are attached to the interior and to every boundary where the
    location applies.

    Args:
        model_part: The model part to convert.
        field_specs: iterable of (variable, data_location) pairs (see
            BuildMesh). Element/condition locations are attached per mesh
            according to each mesh's own source container.
        boundary_sub_model_part_names: names of sub-model-parts to expose as
            DomainMesh boundaries.
        source_container: source container of the interior mesh.
        tessellation_mode: see BuildProvenance. "smallest_id_diagonal" makes
            the boundary meshes' triangulations match the interior mesh's
            boundary faces exactly.
        higher_order_mode: see BuildProvenance.

    Returns:
        (domain_mesh, provenance_maps): provenance_maps maps "interior" and
        each boundary name to its MeshProvenanceMap.
    """
    physicsnemo = _TryImportPhysicsNemo()

    interior_provenance = BuildProvenance(model_part, source_container, tessellation_mode,
                                          higher_order_mode, curved_refinement_levels)
    point_data, cell_data = _CollectFieldData(
        model_part, [(v, loc) for v, loc in field_specs], interior_provenance)
    interior = _MeshFromProvenance(physicsnemo, interior_provenance, point_data, cell_data)
    provenance_maps = {"interior": interior_provenance}

    boundaries = {}
    for name in boundary_sub_model_part_names:
        sub_model_part = model_part.GetSubModelPart(name)
        if sub_model_part.NumberOfConditions() > 0:
            boundary_container = "Conditions"
        elif sub_model_part.NumberOfElements() > 0:
            boundary_container = "Elements"
        else:
            Kratos.Logger.PrintWarning(
                "BuildDomainMesh",
                f"Sub-model-part \"{name}\" has neither conditions nor elements; skipping boundary.")
            continue
        boundary_provenance = BuildProvenance(
            sub_model_part, boundary_container, tessellation_mode,
            higher_order_mode, curved_refinement_levels)
        # Only nodal fields transfer to boundaries generically; entity fields
        # of the interior container do not exist on the boundary container.
        nodal_specs = [(v, loc) for v, loc in field_specs
                       if loc in ("node_historical", "node_non_historical")]
        boundary_point_data, _ = _CollectFieldData(sub_model_part, nodal_specs, boundary_provenance)
        boundaries[name] = _MeshFromProvenance(physicsnemo, boundary_provenance, boundary_point_data, {})
        provenance_maps[name] = boundary_provenance

    domain_mesh = physicsnemo.mesh.DomainMesh(interior=interior, boundaries=boundaries)
    return domain_mesh, provenance_maps


def SaveMesh(mesh, prefix: str):
    """Saves a physicsnemo Mesh/DomainMesh with its native memory-mapped
    on-disk format (.pmsh/.pdmsh directory layout)."""
    return mesh.save(prefix=str(prefix))


def LoadMesh(prefix: str, device=None, domain: bool = False):
    """Loads a mesh saved by SaveMesh (physicsnemo's native format).

    Args:
        prefix: The prefix SaveMesh was given.
        device: Where to place the loaded tensors.
        domain: Load a DomainMesh (interior plus named boundaries) rather
            than a plain Mesh. The two formats are distinct and each
            loader rejects the other's files.
    """
    physicsnemo = _TryImportPhysicsNemo()
    mesh_class = physicsnemo.mesh.DomainMesh if domain else physicsnemo.mesh.Mesh
    if device is None:
        return mesh_class.load(str(prefix))
    return mesh_class.load(str(prefix), device=device)


def _TryImportZarrIo():
    try:
        from physicsnemo.mesh.io import from_zarr, to_zarr
        return to_zarr, from_zarr
    except ImportError as e:
        raise ImportError(
            "Zarr mesh IO requires physicsnemo >= 2.2, which could not be imported. "
            "Install it with e.g. 'pip install -U nvidia-physicsnemo'.") from e


def SaveMeshZarr(mesh, store, chunk_rows: int = 200000, zstd_level: int = 3) -> str:
    """Writes a Mesh or DomainMesh as an AI-ready Zarr store.

    physicsnemo 2.2 gained this, which matters because the alternative -
    physicsnemo-curator's Zarr sink - is a git-only package that downloads
    a Rust toolchain at build time and whose sinks are Mesh-typed, so a
    DomainMesh fails in them. This needs neither.
    """
    to_zarr, _ = _TryImportZarrIo()
    to_zarr(mesh, str(store), chunk_rows=chunk_rows, zstd_level=zstd_level)
    return str(store)


def LoadMeshZarr(store, device=None):
    """Reads a Zarr store written by SaveMeshZarr.

    The store records which of the two it holds, so a DomainMesh comes back
    as a DomainMesh without being asked.
    """
    _, from_zarr = _TryImportZarrIo()
    if device is None:
        return from_zarr(str(store))
    return from_zarr(str(store), device=device)


def ScatterFieldBack(provenance: MeshProvenanceMap,
                     predicted,
                     model_part: Kratos.ModelPart,
                     variable,
                     data_location: str,
                     reduction: str = "mean") -> None:
    """Writes a prediction on the tessellated mesh back onto Kratos entities.

    Args:
        provenance: The MeshProvenanceMap of the tessellation.
        predicted: (P, ...) array for nodal locations (per simplex point) or
            (C, ...) array for element/condition locations (per simplex cell).
        model_part: The model part to write into.
        variable: The Kratos variable to write.
        data_location: "node_historical", "node_non_historical", "element" or
            "condition". Gauss-point write-back is not possible in Kratos
            (GaussPointVariableTensorAdaptor is read-only): write to an
            element variable instead.
        reduction: Reduction for cell fields ("mean", "weighted_mean" using
            sub-cell measures, or "first").
    """
    predicted = numpy.asarray(predicted)

    if data_location in ("node_historical", "node_non_historical"):
        if predicted.shape[0] != provenance.number_of_points:
            raise ValueError(
                f"Nodal prediction has {predicted.shape[0]} rows but the tessellation "
                f"has {provenance.number_of_points} points.")
        node_ids = [node.Id for node in model_part.Nodes]
        nodal_values = provenance.ScatterNodalField(node_ids, predicted)
        tensor_adaptor = GetTensorAdaptor(model_part, data_location, variable, collect=True)
        tensor_adaptor.data[:] = nodal_values.astype(tensor_adaptor.data.dtype, copy=False)
        tensor_adaptor.StoreData()
        return

    if data_location in ("element", "condition"):
        if predicted.shape[0] != provenance.number_of_cells:
            raise ValueError(
                f"Cell prediction has {predicted.shape[0]} rows but the tessellation "
                f"has {provenance.number_of_cells} cells.")
        weights = provenance.ComputeSimplexMeasures() if reduction == "weighted_mean" else None
        entity_ids, values = provenance.AggregateCellField(predicted, reduction, weights)

        tensor_adaptor = GetTensorAdaptor(model_part, data_location, variable, collect=True)
        container = _GetContainer(model_part, provenance.source_container)
        container_ids = numpy.fromiter(
            (entity.Id for entity in container), dtype=numpy.int64, count=len(container))
        if len(container_ids) > 0 and numpy.all(numpy.diff(container_ids) > 0):
            # id-sorted container (the Kratos default): vectorized assignment
            positions = numpy.searchsorted(container_ids, entity_ids)
            valid = (positions < len(container_ids)) & (container_ids[
                numpy.minimum(positions, len(container_ids) - 1)] == entity_ids)
            tensor_adaptor.data[positions[valid]] = values[valid]
        else:  # unsorted container: per-entity fallback
            entity_row = {int(entity_id): row for row, entity_id in enumerate(entity_ids)}
            for row, entity in enumerate(container):
                index = entity_row.get(entity.Id)
                if index is not None:
                    tensor_adaptor.data[row] = values[index]
        tensor_adaptor.StoreData()
        return

    raise ValueError(
        f"Unsupported data location \"{data_location}\" for scatter-back. Gauss-point "
        "write-back is not possible; use \"element\" or \"condition\" instead.")
