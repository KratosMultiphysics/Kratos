"""Particle-graph construction for Lagrangian (Learning-to-Simulate) surrogates.

Particle methods (MPM, SPH, DEM, PFEM) have no persistent element-edge
graph - connectivity is proximity, rebuilt every step. BuildParticleGraph
connects a model part's nodes by radius or kNN and returns exactly
graph_bridge.BuildGraph's contract, so ScatterNodeFeatures, ToPyGGraph and
every downstream idiom work unchanged:

    (node_features (N, F) float64,
     edge_index    (2, E) int64, bidirectional,
     edge_features (E, 4) float64: relative position + distance,
     node_ids      (N,)  int64)

Neighbor search runs through physicsnemo.nn.functional (warp-backed
radius_search/knn) when available and falls back to an exact numpy
brute-force path otherwise ("backend": "numpy" forces it - also the
reference the accelerated path is tested against). "backend": "cupy" runs
that same exact brute-force search on the GPU, needing no physicsnemo.

"box_size" in the connectivity block makes the search PERIODIC (minimum-image
convention): a pair straddling the box boundary is a neighbour with the short
edge vector, which is what a molecular-dynamics cloud (the Lennard-Jones
recipe) needs. Positions may be unwrapped - the images are taken modulo the
box. The periodic search runs on scipy's cKDTree (exact, O(N log N)) for
"auto", and on the brute-force distance matrix for "numpy"/"cupy"; warp's
radius_search is not periodic and is never used with a box.

BuildKinematicFeatures assembles the standard Learning-to-Simulate node
features: the last K velocity states from the historical buffer, oldest
first (matching the TimeSeriesInferenceProcess window convention).

torch/physicsnemo stay optional; the numpy path needs neither.
"""

import numpy

import KratosMultiphysics as Kratos
from KratosMultiphysics.PhysicsNeMoApplication.utilities import array_backend_utils
from KratosMultiphysics.PhysicsNeMoApplication.utilities.tensor_adaptor_dataset_utils import GetTensorAdaptor

_CONNECTIVITY_TYPES = ("radius", "knn")
_BACKENDS = ("auto", "numpy", "cupy")


def _TryImportNeighborSearch():
    try:
        from physicsnemo.nn.functional import knn, radius_search
        return knn, radius_search
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.particle_bridge's accelerated neighbor search requires "
            "physicsnemo (and torch/warp), which could not be imported. Install it with e.g. "
            "'pip install nvidia-physicsnemo', or use the \"numpy\" backend.") from e


def _ReadBox(settings: Kratos.Parameters):
    """The periodic box as a (3,) array, or None (the default: open space).

    "box_size" is [] (non-periodic), [L] (cubic) or [Lx, Ly, Lz].
    """
    box = numpy.asarray(settings["box_size"].GetVector(), dtype=numpy.float64).reshape(-1)
    if box.size == 0:
        return None
    if box.size == 1:
        box = numpy.repeat(box, 3)
    if box.size != 3 or numpy.any(box <= 0.0):
        raise ValueError(
            f"\"box_size\" must be [] (non-periodic), [L] or [Lx, Ly, Lz] with positive "
            f"lengths; got {settings['box_size'].GetVector()}.")
    return box


def _ReadConnectivity(settings: Kratos.Parameters):
    defaults = Kratos.Parameters("""{
        "type"          : "radius",
        "radius"        : 0.015,
        "max_neighbors" : 16,
        "backend"       : "auto",
        "box_size"      : []
    }""")
    settings = settings.Clone()
    settings.ValidateAndAssignDefaults(defaults)
    connectivity_type = settings["type"].GetString()
    if connectivity_type not in _CONNECTIVITY_TYPES:
        raise ValueError(
            f"Unsupported connectivity type \"{connectivity_type}\". "
            f"Use one of {_CONNECTIVITY_TYPES}.")
    backend = settings["backend"].GetString()
    if backend not in _BACKENDS:
        raise ValueError(f"Unsupported backend \"{backend}\". Use one of {_BACKENDS}.")
    radius = settings["radius"].GetDouble()
    if connectivity_type == "radius" and radius <= 0.0:
        raise ValueError(f"\"radius\" must be > 0 [ radius = {radius} ].")
    max_neighbors = settings["max_neighbors"].GetInt()
    if connectivity_type == "knn" and max_neighbors < 1:
        raise ValueError(f"\"max_neighbors\" must be >= 1 [ max_neighbors = {max_neighbors} ].")
    box = _ReadBox(settings)
    if box is not None and connectivity_type == "radius" and radius > 0.5 * box.min():
        raise ValueError(
            f"\"radius\" ({radius}) exceeds half the box ({0.5 * box.min()}); the minimum "
            "image of a pair is then ambiguous.")
    return connectivity_type, radius, max_neighbors, backend, box


def _PeriodicPairs(positions, connectivity_type, radius, max_neighbors, box):
    """Exact periodic neighbour pairs -> (2, E) directed, via scipy's
    cKDTree with boxsize (the tree wants positions inside the box)."""
    from scipy.spatial import cKDTree

    wrapped = numpy.mod(positions, box)
    wrapped[wrapped >= box] -= box[numpy.nonzero(wrapped >= box)[1]]  # exact-boundary guard
    tree = cKDTree(wrapped, boxsize=box)
    n = positions.shape[0]
    if connectivity_type == "radius":
        pairs = tree.query_pairs(radius, output_type="ndarray")  # (P, 2), i < j
        if pairs.size == 0:
            return numpy.zeros((2, 0), dtype=numpy.int64)
        return numpy.concatenate([pairs.T, pairs.T[::-1]], axis=1).astype(numpy.int64)
    k = min(max_neighbors + 1, n)  # +1: the nearest match is self
    _, indices = tree.query(wrapped, k=k)
    indices = numpy.asarray(indices).reshape(n, -1)
    receiver = numpy.repeat(numpy.arange(n), indices.shape[1])
    sender = indices.ravel()
    keep = receiver != sender
    return numpy.stack([sender[keep], receiver[keep]]).astype(numpy.int64)


def _BruteForcePairs(positions, connectivity_type, radius, max_neighbors, backend="numpy",
                     box=None):
    """Exact O(N^2) neighbor pairs -> (2, E) directed [sender, receiver].

    The quadratic distance matrix is the one place in this application where
    a GPU array library pays for itself unambiguously (measured ~3.2x at
    N=2000 and ~3.4x at N=8000, transfers included), because the work grows
    as N^2 while the transfer only grows as N. Selecting "cupy" runs exactly
    this algorithm on the device; the result is returned on the host.
    """
    xp, _ = array_backend_utils.ResolveArrayModule(
        backend, size_hint=positions.shape[0] ** 2)
    positions = xp.asarray(positions)
    n = positions.shape[0]
    deltas = positions[None, :, :] - positions[:, None, :]  # [receiver, sender]
    if box is not None:  # minimum image
        box = xp.asarray(box)
        deltas = deltas - box * xp.round(deltas / box)
    distances = xp.linalg.norm(deltas, axis=-1)
    xp.fill_diagonal(distances, xp.inf)  # no self-edges
    if connectivity_type == "radius":
        receiver, sender = xp.nonzero(distances <= radius)
    else:  # knn: max_neighbors nearest per receiver, then symmetrized
        k = min(max_neighbors, n - 1)
        sender = xp.argsort(distances, axis=1)[:, :k].ravel()
        receiver = xp.repeat(xp.arange(n), k)
    return array_backend_utils.ToHost(xp.stack([sender, receiver])).astype(numpy.int64)


def _AcceleratedPairs(positions, connectivity_type, radius, max_neighbors):
    knn, radius_search = _TryImportNeighborSearch()
    import torch

    points = torch.as_tensor(positions, dtype=torch.float32)
    if connectivity_type == "radius":
        # max_points=None: exact (2, P) [query(receiver), point(sender)] pairs
        pairs = radius_search(points, points, float(radius)).cpu().numpy()
        receiver, sender = pairs[0], pairs[1]
        keep = receiver != sender  # drop self-matches
        return numpy.stack([sender[keep], receiver[keep]]).astype(numpy.int64)
    k = min(max_neighbors + 1, positions.shape[0])  # +1: nearest match is self
    indices, _ = knn(points, points, k)
    indices = indices.cpu().numpy()
    receiver = numpy.repeat(numpy.arange(positions.shape[0]), k)
    sender = indices.ravel()
    keep = receiver != sender
    return numpy.stack([sender[keep], receiver[keep]]).astype(numpy.int64)


def BuildParticleGraphFromPositions(positions, connectivity: Kratos.Parameters):
    """(N, 3) positions -> (edge_index (2, E) bidirectional, edge_features (E, 4)).

    Edge features follow graph_bridge's convention: relative position
    (receiver - sender) plus the distance.
    """
    positions = numpy.asarray(positions, dtype=numpy.float64).reshape(-1, 3)
    connectivity_type, radius, max_neighbors, backend, box = _ReadConnectivity(connectivity)

    if backend in ("numpy", "cupy"):
        # Both force the exact brute-force path; they differ only in where
        # the distance matrix is formed.
        pairs = _BruteForcePairs(positions, connectivity_type, radius, max_neighbors, backend,
                                 box=box)
    elif box is not None:  # periodic: warp's search is not, the KD-tree is
        pairs = _PeriodicPairs(positions, connectivity_type, radius, max_neighbors, box)
    else:  # "auto": the warp-backed neighbour search, exact numpy otherwise
        try:
            pairs = _AcceleratedPairs(positions, connectivity_type, radius, max_neighbors)
        except ImportError:
            pairs = _BruteForcePairs(positions, connectivity_type, radius, max_neighbors)

    # bidirectional with unique edges (radius pairs are symmetric already;
    # knn is not - symmetrize either way for a single canonical result)
    directed = numpy.concatenate([pairs, pairs[::-1]], axis=1)
    directed = numpy.unique(directed, axis=1)
    sender, receiver = directed[0], directed[1]

    relative = positions[receiver] - positions[sender]
    if box is not None:  # the short image, so a boundary-straddling edge is short
        relative = relative - box * numpy.round(relative / box)
    distance = numpy.linalg.norm(relative, axis=1, keepdims=True)
    edge_features = numpy.concatenate([relative, distance], axis=1)
    return directed, edge_features


def BuildParticleGraph(model_part: Kratos.ModelPart, connectivity: Kratos.Parameters,
                       field_specs=()):
    """Builds the proximity graph of a model part's nodes (current positions).

    Returns (node_features, edge_index, edge_features, node_ids) - exactly
    graph_bridge.BuildGraph's contract.
    """
    position_ta = Kratos.TensorAdaptors.NodePositionTensorAdaptor(
        model_part.Nodes, Kratos.Configuration.Current)
    position_ta.CollectData()
    positions = numpy.array(position_ta.data, dtype=numpy.float64).reshape(-1, 3)

    edge_index, edge_features = BuildParticleGraphFromPositions(positions, connectivity)

    node_ids = numpy.fromiter((node.Id for node in model_part.Nodes),
                              dtype=numpy.int64, count=model_part.NumberOfNodes())
    features = []
    for variable_name, data_location in field_specs:
        variable = Kratos.KratosGlobals.GetVariable(variable_name)
        tensor_adaptor = GetTensorAdaptor(model_part, data_location, variable)
        data = numpy.array(tensor_adaptor.data, dtype=numpy.float64)
        features.append(data.reshape(data.shape[0], -1))
    node_features = (numpy.concatenate(features, axis=1) if features
                     else numpy.zeros((len(node_ids), 0)))
    return node_features, edge_index, edge_features, node_ids


def BuildKinematicFeatures(model_part: Kratos.ModelPart, history_size: int,
                           velocity_variable=None):
    """(N, K*3) velocity history from the historical buffer, oldest first.

    The model part's buffer size must be >= history_size (the standard
    Learning-to-Simulate node features; append node-type one-hots yourself).
    """
    if velocity_variable is None:
        velocity_variable = Kratos.VELOCITY
    if history_size < 1:
        raise ValueError(f"history_size must be >= 1 [ history_size = {history_size} ].")
    if model_part.GetBufferSize() < history_size:
        raise ValueError(
            f"The model part's buffer size ({model_part.GetBufferSize()}) is smaller than "
            f"history_size ({history_size}); increase the buffer size.")

    n = model_part.NumberOfNodes()
    blocks = []
    for age in range(history_size - 1, -1, -1):  # oldest first
        block = numpy.empty((n, 3))
        for row, node in enumerate(model_part.Nodes):
            block[row] = node.GetSolutionStepValue(velocity_variable, age)
        blocks.append(block)
    return numpy.concatenate(blocks, axis=1)
