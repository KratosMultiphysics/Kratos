"""Bridge between Kratos meshes and graph representations for GNNs.

Graph neural networks (e.g. physicsnemo.models.meshgraphnet.MeshGraphNet)
operate directly on the mesh graph — the native representation of a finite
element model, with no voxelization loss. This module extracts
(node_features, edge_index, edge_features) from a ModelPart and writes
per-node predictions back.

Edge convention: unique geometric element edges (corner-node pairs along each
geometry's edges — no cell diagonals), both directions included (the
MeshGraphNet bidirectional convention). Edge features are the standard
MeshGraphNet encoding: relative position (x_j - x_i, 3 components) plus the
Euclidean distance (1 component).

Pure Kratos + numpy at module scope; torch_geometric is imported lazily in
ToPyGGraph only.
"""

import numpy

import KratosMultiphysics as Kratos
from KratosMultiphysics.PhysicsNeMoApplication.utilities import array_backend_utils
from KratosMultiphysics.PhysicsNeMoApplication.utilities.tensor_adaptor_dataset_utils import (
    GetTensorAdaptor, RowsOfIds)

_GEOMETRY_TYPE = Kratos.GeometryData.KratosGeometryType

# Corner-edge tables in local corner indices (higher-order types share the
# table of their linear counterpart; corners come first in Kratos ordering).
_EDGE_TABLES = {
    # triangles
    _GEOMETRY_TYPE.Kratos_Triangle2D3: ((0, 1), (1, 2), (2, 0)),
    _GEOMETRY_TYPE.Kratos_Triangle3D3: ((0, 1), (1, 2), (2, 0)),
    _GEOMETRY_TYPE.Kratos_Triangle2D6: ((0, 1), (1, 2), (2, 0)),
    _GEOMETRY_TYPE.Kratos_Triangle3D6: ((0, 1), (1, 2), (2, 0)),
    # quadrilaterals
    _GEOMETRY_TYPE.Kratos_Quadrilateral2D4: ((0, 1), (1, 2), (2, 3), (3, 0)),
    _GEOMETRY_TYPE.Kratos_Quadrilateral3D4: ((0, 1), (1, 2), (2, 3), (3, 0)),
    _GEOMETRY_TYPE.Kratos_Quadrilateral2D8: ((0, 1), (1, 2), (2, 3), (3, 0)),
    _GEOMETRY_TYPE.Kratos_Quadrilateral2D9: ((0, 1), (1, 2), (2, 3), (3, 0)),
    _GEOMETRY_TYPE.Kratos_Quadrilateral3D8: ((0, 1), (1, 2), (2, 3), (3, 0)),
    _GEOMETRY_TYPE.Kratos_Quadrilateral3D9: ((0, 1), (1, 2), (2, 3), (3, 0)),
    # tetrahedra
    _GEOMETRY_TYPE.Kratos_Tetrahedra3D4: ((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)),
    _GEOMETRY_TYPE.Kratos_Tetrahedra3D10: ((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)),
    # hexahedra (bottom ring, top ring, verticals)
    _GEOMETRY_TYPE.Kratos_Hexahedra3D8: (
        (0, 1), (1, 2), (2, 3), (3, 0),
        (4, 5), (5, 6), (6, 7), (7, 4),
        (0, 4), (1, 5), (2, 6), (3, 7)),
    _GEOMETRY_TYPE.Kratos_Hexahedra3D20: (
        (0, 1), (1, 2), (2, 3), (3, 0),
        (4, 5), (5, 6), (6, 7), (7, 4),
        (0, 4), (1, 5), (2, 6), (3, 7)),
    _GEOMETRY_TYPE.Kratos_Hexahedra3D27: (
        (0, 1), (1, 2), (2, 3), (3, 0),
        (4, 5), (5, 6), (6, 7), (7, 4),
        (0, 4), (1, 5), (2, 6), (3, 7)),
    # prisms/wedges (bottom triangle, top triangle, verticals)
    _GEOMETRY_TYPE.Kratos_Prism3D6: (
        (0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4), (2, 5)),
    _GEOMETRY_TYPE.Kratos_Prism3D15: (
        (0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4), (2, 5)),
    # pyramids (base ring, apex edges)
    _GEOMETRY_TYPE.Kratos_Pyramid3D5: (
        (0, 1), (1, 2), (2, 3), (3, 0), (0, 4), (1, 4), (2, 4), (3, 4)),
    _GEOMETRY_TYPE.Kratos_Pyramid3D13: (
        (0, 1), (1, 2), (2, 3), (3, 0), (0, 4), (1, 4), (2, 4), (3, 4)),
    # lines
    _GEOMETRY_TYPE.Kratos_Line2D2: ((0, 1),),
    _GEOMETRY_TYPE.Kratos_Line3D2: ((0, 1),),
}

_NODAL_LOCATIONS = ("node_historical", "node_non_historical")


def _TryImportPyG():
    try:
        import torch_geometric.data
        return torch_geometric
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.graph_bridge.ToPyGGraph requires torch_geometric, which "
            "could not be imported. Install it with e.g. 'pip install torch_geometric'.") from e


def _GetContainer(model_part: Kratos.ModelPart, source_container: str):
    if source_container == "Elements":
        return model_part.Elements
    if source_container == "Conditions":
        return model_part.Conditions
    raise ValueError(f"Unsupported source container \"{source_container}\". Use \"Elements\" or \"Conditions\".")


def BuildGraph(model_part: Kratos.ModelPart, field_specs=(), source_container: str = "Elements"):
    """Extracts the mesh graph of a model part.

    Args:
        model_part: The model part.
        field_specs: iterable of (variable_name, data_location) pairs (nodal
            locations only) forming the node features.
        source_container: "Elements" (default) or "Conditions".

    Returns:
        (node_features, edge_index, edge_features, node_ids):
        node_features (N, F) float64 (F = 0 for an empty spec);
        edge_index (2, E) int64 rows into the node arrays, bidirectional;
        edge_features (E, 4) float64: relative position + distance;
        node_ids (N,) int64 Kratos node ids per row.
    """
    node_ids = numpy.fromiter((node.Id for node in model_part.Nodes), dtype=numpy.int64)
    node_row = {int(node_id): row for row, node_id in enumerate(node_ids)}

    undirected = set()
    for entity in _GetContainer(model_part, source_container):
        geometry = entity.GetGeometry()
        geometry_type = geometry.GetGeometryType()
        if geometry_type not in _EDGE_TABLES:
            raise RuntimeError(
                f"Unsupported geometry type for graph extraction: {geometry_type}.")
        entity_node_ids = [node.Id for node in geometry]
        for a, b in _EDGE_TABLES[geometry_type]:
            i, j = node_row[entity_node_ids[a]], node_row[entity_node_ids[b]]
            undirected.add((i, j) if i < j else (j, i))

    pairs = numpy.array(sorted(undirected), dtype=numpy.int64).reshape(-1, 2)
    edge_index = numpy.concatenate([pairs.T, pairs.T[::-1]], axis=1)  # (2, 2*|undirected|)

    edge_features = ComputeEdgeFeatures(model_part, edge_index)
    node_features = GatherNodeFeatures(model_part, field_specs, len(node_ids))

    return node_features, edge_index, edge_features, node_ids


def ComputeEdgeFeatures(model_part: Kratos.ModelPart, edge_index, backend="numpy"):
    """Relative position + distance for a *given* edge index.

    Split out of `BuildGraph` because the edge index is topology and the
    features are geometry: a deforming mesh needs these recomputed every
    step, but re-extracting the edge set to get them costs orders of
    magnitude more than the arithmetic (see benchmarks/benchmark_bridges.py
    and GraphInferenceProcess, which caches the index).

    Args:
        model_part: The model part the edge index indexes into.
        edge_index: (2, E) int64 rows into the node arrays.
        backend: "numpy" (default), "cupy" or "auto". The CuPy path measured
            ~2.7x at 200k edges and ~3.5x at 1.2M, transfers included; below
            the size threshold it falls back to numpy on its own.

    Returns:
        (E, 4) float64: relative position + distance.
    """
    position_ta = Kratos.TensorAdaptors.NodePositionTensorAdaptor(model_part.Nodes, Kratos.Configuration.Current)
    position_ta.CollectData()
    xp, _ = array_backend_utils.ResolveArrayModule(
        backend, size_hint=numpy.asarray(edge_index).shape[-1])
    coordinates = xp.asarray(position_ta.data)
    edge_index = xp.asarray(edge_index)
    relative = coordinates[edge_index[1]] - coordinates[edge_index[0]]
    distance = xp.linalg.norm(relative, axis=1, keepdims=True)
    return array_backend_utils.ToHost(xp.concatenate([relative, distance], axis=1))


def GatherNodeFeatures(model_part: Kratos.ModelPart, field_specs=(), num_nodes=None):
    """Node features alone, without extracting the graph.

    Rows follow `model_part.Nodes` order, which is the order `BuildGraph`
    reports as `node_ids`, so the two stay aligned.

    Args:
        model_part: The model part.
        field_specs: iterable of (variable_name, data_location) pairs
            (nodal locations only).
        num_nodes: Row count for the empty-spec case; counted if omitted.

    Returns:
        (N, F) float64, F = 0 for an empty spec.
    """
    features = []
    for variable_name, data_location in field_specs:
        if data_location not in _NODAL_LOCATIONS:
            raise ValueError(
                f"Graph node features support nodal locations only "
                f"({', '.join(_NODAL_LOCATIONS)}), got \"{data_location}\".")
        variable = Kratos.KratosGlobals.GetVariable(variable_name)
        tensor_adaptor = GetTensorAdaptor(model_part, data_location, variable)
        data = numpy.array(tensor_adaptor.data)
        features.append(data.reshape(data.shape[0], -1))
    if features:
        return numpy.concatenate(features, axis=1)
    if num_nodes is None:
        num_nodes = model_part.NumberOfNodes()
    return numpy.zeros((num_nodes, 0))


def BuildScatterRows(model_part: Kratos.ModelPart, node_ids):
    """node_ids -> row indices into the model part's nodal arrays.

    Topology, so a caller stepping a static mesh can build this once and
    hand it to `ScatterNodeFeatures` instead of paying the O(N) dict per
    step.
    """
    part_ids = numpy.fromiter((node.Id for node in model_part.Nodes),
                              dtype=numpy.int64, count=model_part.NumberOfNodes())
    # the shared searchsorted lookup (utilities.tensor_adaptor_dataset_utils
    # .RowsOfIds) rather than a {id: row} dict plus a generator; this is the
    # fallback path when a caller has not cached the mapping
    return RowsOfIds(part_ids, numpy.asarray(node_ids, dtype=numpy.int64).ravel())


def ScatterNodeFeatures(model_part: Kratos.ModelPart, node_ids, values,
                        output_field_specs, rows=None) -> None:
    """Writes per-node values (aligned with node_ids) back onto the model part.

    Channels are split across the output fields by each field's per-node
    width, exactly like InferenceProcess.

    Args:
        rows: Optional precomputed `BuildScatterRows` result. The mapping is
            topological, so a caller that knows the mesh has not changed can
            reuse it across steps.
    """
    values = numpy.asarray(values)
    if rows is None:
        rows = BuildScatterRows(model_part, node_ids)

    offset = 0
    for variable_name, data_location in output_field_specs:
        if data_location not in _NODAL_LOCATIONS:
            raise ValueError(
                f"Graph scatter supports nodal locations only ({', '.join(_NODAL_LOCATIONS)}), "
                f"got \"{data_location}\".")
        variable = Kratos.KratosGlobals.GetVariable(variable_name)
        tensor_adaptor = GetTensorAdaptor(model_part, data_location, variable)
        width = int(numpy.prod(tensor_adaptor.data.shape[1:], dtype=int))
        chunk = numpy.zeros(tensor_adaptor.data.shape)
        chunk.reshape(chunk.shape[0], -1)[rows] = values[:, offset:offset + width]
        tensor_adaptor.data[:] = chunk
        tensor_adaptor.StoreData()
        offset += width
    if offset != values.shape[1]:
        raise ValueError(
            f"Prediction has {values.shape[1]} channels but the output fields consume {offset}.")


def ToPyGGraph(edge_index, num_nodes: int, positions=None):
    """Builds the torch_geometric graph object MeshGraphNet's forward expects.

    Args:
        edge_index: (2, E) int64 array (as returned by BuildGraph).
        num_nodes: Number of graph nodes.
        positions: Optional (N, D) node coordinates stored as `Data.pos`.
            BiStrideMeshGraphNet reads `graph.pos` unconditionally (it has no
            fallback), so the bistride interface requires this; the plain
            MeshGraphNet path ignores it.

    Returns:
        torch_geometric.data.Data with edge_index, num_nodes and optionally pos.
    """
    torch_geometric = _TryImportPyG()
    import torch  # torch_geometric guarantees torch is present
    graph = torch_geometric.data.Data(
        edge_index=torch.from_numpy(numpy.ascontiguousarray(edge_index)),
        num_nodes=num_nodes)
    if positions is not None:
        graph.pos = torch.from_numpy(
            numpy.ascontiguousarray(numpy.asarray(positions, dtype=numpy.float64)))
    return graph


def NodePositions(model_part: Kratos.ModelPart) -> numpy.ndarray:
    """(N, 3) current nodal coordinates in BuildGraph's node-row order."""
    position_ta = Kratos.TensorAdaptors.NodePositionTensorAdaptor(
        model_part.Nodes, Kratos.Configuration.Current)
    position_ta.CollectData()
    return numpy.array(position_ta.data)


def _ParityCoarseSelection(adjacency, positions, component_rows):
    """BSMS node selection inside one connected component.

    Breadth-first search from the node nearest the component's geometric
    centre two-colours it; keeping one BFS parity class halves the node count
    while guaranteeing the kept set is a maximal independent-ish spread. The
    smaller class is kept (upstream's rule), which also keeps the seed's own
    class when the component is a single node.
    """
    from collections import deque

    centre = positions[component_rows].mean(axis=0)
    seed = component_rows[int(numpy.argmin(
        numpy.linalg.norm(positions[component_rows] - centre, axis=1)))]

    depth = {int(seed): 0}
    queue = deque([int(seed)])
    indptr, indices = adjacency.indptr, adjacency.indices
    while queue:
        row = queue.popleft()
        for neighbour in indices[indptr[row]:indptr[row + 1]]:
            neighbour = int(neighbour)
            if neighbour not in depth:
                depth[neighbour] = depth[row] + 1
                queue.append(neighbour)

    even = [row for row, level in depth.items() if level % 2 == 0]
    odd = [row for row, level in depth.items() if level % 2 == 1]
    if not odd:
        return even
    return even if len(even) <= len(odd) else odd


def BuildBistrideHierarchy(edge_index, num_nodes: int, positions, num_levels: int = 1):
    """Multiscale pooling tables for BiStrideMeshGraphNet.

    Implements the BSMS hierarchy (Cao et al., "Efficient Learning of Mesh-Based
    Physical Simulation with Bi-Stride Multi-Scale Graph Neural Network"): each
    level two-colours every connected component by breadth-first distance from
    a geometric-centre seed, keeps one parity class, and connects the survivors
    through the squared adjacency (two-hop neighbours in the parent level).

    physicsnemo ships the same algorithm in
    physicsnemo.datapipes.gnn.bsms.BistrideMultiLayerGraph, but that path calls
    sparse_dot_mkl for the adjacency square and raises without it (and without
    a working MKL runtime), so this implementation uses scipy - a hard Kratos
    core dependency - instead. It also sets the adjacency diagonal before
    squaring: upstream's message-passing helper derives its node count from
    max(edge index), so a level whose highest-numbered node sources no edge
    would otherwise break, and isolated nodes would divide by a zero degree.

    Args:
        edge_index: (2, E) int64 edge rows of the FULL graph (bidirectional).
        num_nodes: Node count of the full graph.
        positions: (N, D) node coordinates (the seed choice is geometric).
        num_levels: Number of coarse levels to build (the model's
            `num_mesh_levels`).

    Returns:
        (ms_edges, ms_ids): ms_edges holds num_levels + 1 arrays of shape
        (2, E_i) int64, each renumbered into ITS OWN level's node space
        (entry 0 is the full graph); ms_ids holds num_levels arrays of shape
        (N_{i+1},) int64 selecting which rows of level i survive into level
        i + 1. The lengths differ by one - that is the model's contract, not
        an oversight.
    """
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components

    if num_levels < 1:
        raise ValueError(f"num_levels must be >= 1, got {num_levels}.")
    edge_index = numpy.asarray(edge_index, dtype=numpy.int64).reshape(2, -1)
    positions = numpy.asarray(positions, dtype=numpy.float64)

    level_edges = edge_index
    level_positions = positions
    level_nodes = int(num_nodes)

    ms_edges, ms_ids = [], []
    for level in range(num_levels + 1):
        ms_edges.append(numpy.ascontiguousarray(level_edges))
        if level == num_levels:
            break

        data = numpy.ones(level_edges.shape[1], dtype=numpy.int8)
        adjacency = coo_matrix((data, (level_edges[0], level_edges[1])),
                               shape=(level_nodes, level_nodes)).tocsr()
        adjacency.data[:] = 1

        _, labels = connected_components(adjacency, directed=False)
        kept = []
        for label in range(labels.max() + 1 if level_nodes else 0):
            component_rows = numpy.flatnonzero(labels == label)
            kept.extend(_ParityCoarseSelection(adjacency, level_positions, component_rows))
        kept = numpy.sort(numpy.asarray(kept, dtype=numpy.int64))

        if kept.size < 2:
            raise ValueError(
                f"The bistride hierarchy collapsed to {kept.size} node(s) at level "
                f"{level + 1} of {num_levels}; use fewer levels or a finer mesh.")

        # self-loops before squaring: two-hop reachability must include the
        # kept node's own one-hop neighbours that were dropped
        with_loops = adjacency.tolil()
        with_loops.setdiag(1)
        squared = (with_loops.tocsr() @ with_loops.tocsr()).tocoo()

        renumber = numpy.full(level_nodes, -1, dtype=numpy.int64)
        renumber[kept] = numpy.arange(kept.size, dtype=numpy.int64)
        rows, columns = renumber[squared.row], renumber[squared.col]
        keep_edge = (rows >= 0) & (columns >= 0) & (rows != columns)

        ms_ids.append(kept)
        level_edges = numpy.stack([rows[keep_edge], columns[keep_edge]], axis=0)
        level_edges = numpy.unique(level_edges, axis=1)
        level_positions = level_positions[kept]
        level_nodes = int(kept.size)

    return ms_edges, ms_ids


def BuildWorldEdges(model_part: Kratos.ModelPart, connectivity: Kratos.Parameters):
    """Proximity ("world") edges for HybridMeshGraphNet.

    Mesh edges follow the element topology; world edges connect nodes that are
    merely CLOSE in space (contact, free surfaces, self-approach) - the second
    edge type of the MeshGraphNet paper. Built by the shipped particle-bridge
    neighbour search, which already emits this module's edge convention.

    Args:
        model_part: The model part (current coordinates).
        connectivity: particle_bridge connectivity block, e.g.
            {"type": "radius", "radius": 0.1, "backend": "auto"}.

    Returns:
        (edge_index (2, E) int64, edge_features (E, 4) float64).
    """
    from KratosMultiphysics.PhysicsNeMoApplication.bridges import particle_bridge
    return particle_bridge.BuildParticleGraphFromPositions(
        NodePositions(model_part), connectivity)


def ConcatenateEdgeSets(mesh_edge_index, world_edge_index):
    """The combined edge_index HybridMeshGraphNet's single graph object needs.

    The model takes mesh and world edge FEATURES as two separate tensors but
    one graph, and splits the concatenated edge set positionally by row count
    (`efeat[:len(mesh)]` / `efeat[len(mesh):]`). Mesh edges must therefore come
    first - getting it backwards runs silently and trains nonsense, which is
    exactly why this helper exists instead of an inline `numpy.concatenate`.
    """
    mesh_edge_index = numpy.asarray(mesh_edge_index, dtype=numpy.int64).reshape(2, -1)
    world_edge_index = numpy.asarray(world_edge_index, dtype=numpy.int64).reshape(2, -1)
    return numpy.concatenate([mesh_edge_index, world_edge_index], axis=1)
