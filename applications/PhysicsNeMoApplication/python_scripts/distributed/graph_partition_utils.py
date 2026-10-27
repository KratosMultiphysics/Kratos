"""Per-rank mesh subgraphs with halos, for distributed GNN training.

`graph_bridge.BuildGraph` is correct on a serial model part and quietly
wrong on a distributed one. Two things go wrong at a partition interface:

- rank subgraphs OVERLAP - each rank's `model_part.Nodes` is local plus
  ghost, so interface nodes and their edges appear on both sides and the
  per-rank node counts sum to more than the global count;
- neighbourhoods are TRUNCATED - a rank holds only its own elements, so an
  interface node is missing the edges contributed by elements on the far
  side. Measured against a serial reference, a one-hop aggregation is wrong
  at every interface node, *including at nodes the rank owns*.

This module fixes both. Each rank builds a subgraph over the nodes it owns
plus a halo ring of the elements touching them, so:

- owned sets partition the global node set exactly (no double counting), and
- a one-hop aggregation matches the serial reference at every owned node.

An L-layer message-passing network needs L halo rings.

Two rules worth stating because the obvious alternatives fail:

- **Edge ownership follows ELEMENT ownership, never PARTITION_INDEX.**
  Elements are uniquely owned; nodes are not. Deriving edge ownership from
  a `min(PARTITION_INDEX)` rule silently DROPS edges whose endpoints are
  owned by different ranks, because such an edge may not be present on the
  lower-ranked owner at all.
- **Kratos provides no element halo.** `GhostMesh()` carries ghost nodes but
  zero ghost elements, so the neighbouring elements must be exchanged
  explicitly rather than read off the communicator.

Scope: this is DATA parallelism - each rank trains the same model on its own
subgraph and gradients are averaged. Domain-parallel sharding of a single
tensor's activations (physicsnemo's ShardTensor story) lives in
domain_parallel_utils, which builds the device mesh this module does not.

torch is imported lazily; the partitioning itself is pure Kratos + numpy.
"""

import numpy

import KratosMultiphysics as Kratos

from KratosMultiphysics.PhysicsNeMoApplication.bridges import graph_bridge
from KratosMultiphysics.PhysicsNeMoApplication.utilities.tensor_adaptor_dataset_utils import GetTensorAdaptor

_PARTITION_UNSAFE_INTERFACES = ("bistride", "hybrid")


def _TryImportTorch():
    try:
        import torch
        return torch
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.graph_partition_utils requires torch, which could "
            "not be imported. Install it with e.g. 'pip install torch'.") from e


def CheckPartitionSafeInterface(model_part: Kratos.ModelPart, model_interface: str) -> None:
    """Refuses the graph interfaces that cannot be partitioned.

    "bistride" builds its coarse levels from a global connected-components
    pass and a geometric-centre seed, so per-rank hierarchies do not
    correspond to a partition of the global one. "hybrid" adds proximity
    edges by a radius search that ignores the partition entirely, so a
    node's true neighbours may live on another rank and would be silently
    missing. Both would train on quietly wrong graphs.
    """
    if model_interface in _PARTITION_UNSAFE_INTERFACES and model_part.IsDistributed():
        raise ValueError(
            f"The \"{model_interface}\" interface is not partition-safe: its graph is "
            "built from global information (the bistride hierarchy from a global "
            "connected-components pass, hybrid's world edges from a partition-blind "
            "radius search), so a per-rank subgraph would silently differ from the "
            "serial one. Use \"meshgraphnet\"/\"meshgraphkan\" for distributed runs.")


def _LocalElementConnectivity(model_part: Kratos.ModelPart, source_container: str):
    """(entity_id, node_ids) for this rank's OWNED entities."""
    communicator = model_part.GetCommunicator()
    local_mesh = communicator.LocalMesh()
    container = (local_mesh.Elements if source_container == "Elements"
                 else local_mesh.Conditions)
    ids, connectivity = [], []
    for entity in container:
        ids.append(entity.Id)
        connectivity.append([node.Id for node in entity.GetGeometry()])
    return ids, connectivity


def _AllGatherConnectivity(model_part: Kratos.ModelPart, source_container: str,
                           data_communicator):
    """Every rank's owned element connectivity, on every rank.

    Kratos has no element halo (GhostMesh carries zero elements), so the
    neighbouring connectivity has to be exchanged. This gathers all of it,
    which is exact and simple; a production-scale version would exchange
    only with `NeighbourIndices()` colours.
    """
    ids, connectivity = _LocalElementConnectivity(model_part, source_container)
    width = len(connectivity[0]) if connectivity else 0
    widths = data_communicator.AllGatherInts([width])
    width = max(widths) if widths else 0
    if width == 0:
        return [], []

    flat = [int(node_id) for row in connectivity for node_id in row]
    gathered_ids = data_communicator.AllGathervInts([int(i) for i in ids])
    gathered_flat = data_communicator.AllGathervInts(flat)

    all_ids, all_connectivity = [], []
    for rank_ids, rank_flat in zip(gathered_ids, gathered_flat):
        rank_ids = list(rank_ids)
        rank_flat = list(rank_flat)
        for index, entity_id in enumerate(rank_ids):
            all_ids.append(int(entity_id))
            all_connectivity.append(rank_flat[index * width:(index + 1) * width])
    return all_ids, all_connectivity


def BuildHaloSubgraph(model_part: Kratos.ModelPart, num_halo_rings: int = 1,
                      source_container: str = "Elements", field_specs=()):
    """A per-rank subgraph over owned nodes plus `num_halo_rings` of halo.

    Args:
        model_part: A distributed model part (a serial one yields the whole
            graph with every node owned, so callers need no special case).
        num_halo_rings: Halo depth. An L-layer message-passing network needs
            L rings for its owned-node outputs to match a serial run.
        source_container: "Elements" (default) or "Conditions".
        field_specs: iterable of (variable_name, data_location) node features.

    Returns:
        (node_features, edge_index, edge_features, node_ids, owned_mask):
        graph_bridge.BuildGraph's contract plus a boolean (N,) mask marking
        the rows this rank owns. Only owned rows are this rank's to predict
        or scatter back; halo rows exist to make their neighbourhoods whole.
    """
    if num_halo_rings < 1:
        raise ValueError(f"num_halo_rings must be >= 1, got {num_halo_rings}.")

    if not model_part.IsDistributed():
        features, edge_index, edge_features, node_ids = graph_bridge.BuildGraph(
            model_part, field_specs, source_container)
        return (features, edge_index, edge_features, node_ids,
                numpy.ones(len(node_ids), dtype=bool))

    communicator = model_part.GetCommunicator()
    data_communicator = communicator.GetDataCommunicator()
    owned_ids = {node.Id for node in communicator.LocalMesh().Nodes}

    all_ids, all_connectivity = _AllGatherConnectivity(
        model_part, source_container, data_communicator)
    # A halo deeper than one ring reaches nodes outside Kratos's ghost layer,
    # so their coordinates and features must be exchanged too - the ghost
    # mesh alone is not enough.
    global_coordinates, global_features = _AllGatherNodeData(
        model_part, field_specs, data_communicator)

    # grow the node set ring by ring: an element joins when it touches the
    # set, and contributes its nodes to the next ring
    selected_nodes = set(owned_ids)
    selected_elements = set()
    for _ in range(num_halo_rings):
        newly_selected = []
        for index, connectivity in enumerate(all_connectivity):
            if index in selected_elements:
                continue
            if any(node_id in selected_nodes for node_id in connectivity):
                selected_elements.add(index)
                newly_selected.append(connectivity)
        for connectivity in newly_selected:
            selected_nodes.update(connectivity)

    # rows in ascending id order, matching BuildGraph's convention
    node_ids = numpy.array(sorted(selected_nodes), dtype=numpy.int64)
    row_of = {int(node_id): row for row, node_id in enumerate(node_ids)}
    owned_mask = numpy.array([int(node_id) in owned_ids for node_id in node_ids], dtype=bool)

    undirected = set()
    for index in sorted(selected_elements):
        connectivity = all_connectivity[index]
        corners = [node_id for node_id in connectivity if node_id in row_of]
        for first in range(len(corners)):
            for second in range(first + 1, len(corners)):
                a, b = row_of[corners[first]], row_of[corners[second]]
                undirected.add((a, b) if a < b else (b, a))
    pairs = numpy.array(sorted(undirected), dtype=numpy.int64).reshape(-1, 2)
    edge_index = numpy.concatenate([pairs.T, pairs.T[::-1]], axis=1)

    coordinates = numpy.array([global_coordinates[int(i)] for i in node_ids], dtype=float)
    relative = coordinates[edge_index[1]] - coordinates[edge_index[0]]
    distance = numpy.linalg.norm(relative, axis=1, keepdims=True)
    edge_features = numpy.concatenate([relative, distance], axis=1)

    if field_specs:
        node_features = numpy.array(
            [global_features[int(i)] for i in node_ids], dtype=float)
    else:
        node_features = numpy.zeros((len(node_ids), 0))
    return node_features, edge_index, edge_features, node_ids, owned_mask


def _AllGatherNodeData(model_part: Kratos.ModelPart, field_specs, data_communicator):
    """Every OWNED node's coordinates and features, keyed by id, on every rank.

    Owned nodes are gathered (not all local nodes) so each node is
    contributed exactly once. As with the connectivity, this is the exact
    and simple exchange; a production version would trade only with
    neighbouring ranks.
    """
    communicator = model_part.GetCommunicator()
    owned = list(communicator.LocalMesh().Nodes)
    owned_ids = [node.Id for node in owned]
    coordinates = [value for node in owned for value in (node.X, node.Y, node.Z)]

    gathered_ids = data_communicator.AllGathervInts([int(i) for i in owned_ids])
    gathered_coordinates = data_communicator.AllGathervDoubles(
        [float(v) for v in coordinates])

    global_coordinates = {}
    for rank_ids, rank_values in zip(gathered_ids, gathered_coordinates):
        rank_ids = list(rank_ids)
        rank_values = list(rank_values)
        for index, node_id in enumerate(rank_ids):
            global_coordinates[int(node_id)] = rank_values[index * 3:index * 3 + 3]

    global_features = {}
    if field_specs:
        row_of = {node.Id: row for row, node in enumerate(model_part.Nodes)}
        blocks = []
        for variable_name, data_location in field_specs:
            variable = Kratos.KratosGlobals.GetVariable(variable_name)
            tensor_adaptor = GetTensorAdaptor(model_part, data_location, variable)
            data = numpy.array(tensor_adaptor.data)
            blocks.append(data.reshape(data.shape[0], -1))
        stacked = numpy.concatenate(blocks, axis=1)
        width = int(stacked.shape[1])
        owned_rows = [row_of[node_id] for node_id in owned_ids]
        flat = [float(v) for v in stacked[owned_rows].reshape(-1)]
        gathered_features = data_communicator.AllGathervDoubles(flat)
        for rank_ids, rank_values in zip(gathered_ids, gathered_features):
            rank_ids = list(rank_ids)
            rank_values = list(rank_values)
            for index, node_id in enumerate(rank_ids):
                global_features[int(node_id)] = rank_values[index * width:(index + 1) * width]
    return global_coordinates, global_features


def _TryImportAutogradCollectives():
    try:
        from physicsnemo.distributed.autograd import indexed_all_to_all_v
        return indexed_all_to_all_v
    except ImportError as e:
        raise ImportError(
            "The differentiable halo exchange requires physicsnemo, which could not be "
            "imported. Install it with e.g. 'pip install -U nvidia-physicsnemo'.") from e


class HaloExchangePlan:
    """Who sends which OWNED rows to whom, so halo features can travel
    through torch autograd instead of being copied in up front.

    BuildHaloSubgraph exchanges node features eagerly over the Kratos
    DataCommunicator: correct for a forward pass, but a copy with no
    gradient path. A model whose INPUT features are themselves being
    learned - an encoder upstream of the graph, or features that are
    optimization variables - needs the gradient a halo copy receives to flow
    back to the rank that owns the original. This plan drives
    physicsnemo.distributed.autograd.indexed_all_to_all_v, whose backward
    does exactly that: it SUMS the gradient of every copy onto the owner.

    The topology half stays eager (connectivity carries no gradient); only
    the features move through torch.

    PADDING, and why it is not optional. A halo exchange is ragged by
    nature: a rank sends many rows to the neighbour it shares an interface
    with and none to a rank on the far side of the mesh. Gloo's all-to-all
    requires every block in the list to have the SAME shape and rejects
    ragged ones outright ("invalid tensor size at index 1"), so every block
    is padded to one common size and sliced back on receipt. The padding
    rows are repeats of an existing index, and since the received copies of
    them are discarded, their gradient is exactly zero - the property the
    exchange exists for is untouched. A NCCL transport would not need this,
    but NCCL needs one GPU per rank, which is why the padded path is the
    one that is actually tested here.

    Attributes:
        send_rows: per destination rank, a LongTensor of rows into this
            rank's owned features, in ascending node-id order, padded to
            block_size.
        true_sizes: true_sizes[i][j] = how many rows rank i really sends to
            rank j, before padding.
        block_size: the padded block size every pair exchanges.
        owned_rows / halo_rows: where owned and received rows land in the
            subgraph's node_ids order.
    """

    def __init__(self, node_ids, owned_mask, data_communicator) -> None:
        torch = _TryImportTorch()
        node_ids = numpy.asarray(node_ids, dtype=numpy.int64)
        owned_mask = numpy.asarray(owned_mask, dtype=bool)
        self.size = data_communicator.Size()
        self.rank = data_communicator.Rank()
        self.n_nodes = len(node_ids)
        self.n_owned = int(owned_mask.sum())

        owned_ids = numpy.sort(node_ids[owned_mask])
        needed_ids = numpy.sort(node_ids[~owned_mask])
        gathered_owned = [numpy.asarray(list(ids), dtype=numpy.int64) for ids in
                          data_communicator.AllGathervInts([int(i) for i in owned_ids])]
        gathered_needed = [numpy.asarray(list(ids), dtype=numpy.int64) for ids in
                           data_communicator.AllGathervInts([int(i) for i in needed_ids])]

        # true_sizes[i][j]: what rank i owns that rank j's halo needs
        self.true_sizes = [
            [int(numpy.intersect1d(gathered_owned[i], gathered_needed[j]).size)
             for j in range(self.size)]
            for i in range(self.size)]
        self.block_size = max(
            (count for row in self.true_sizes for count in row), default=0)
        # what upstream is told: one uniform block per pair (see PADDING)
        self.sizes = [[self.block_size] * self.size for _ in range(self.size)]

        if self.block_size and self.n_owned == 0:
            raise ValueError(
                "This rank owns no node but the exchange needs padding rows taken from "
                "its own features. A partition that leaves a rank empty cannot use the "
                "differentiable halo exchange.")

        self.send_rows = []
        for destination in range(self.size):
            shared = numpy.intersect1d(owned_ids, gathered_needed[destination])
            rows = numpy.searchsorted(owned_ids, shared).astype(numpy.int64)
            if self.block_size:
                # pad with an existing row: its received copies are dropped,
                # so it contributes exactly zero gradient
                padding = numpy.zeros(self.block_size - rows.size, dtype=numpy.int64)
                rows = numpy.concatenate([rows, padding])
            self.send_rows.append(torch.from_numpy(rows))

        # what arrives here is one padded block per source, each holding its
        # true rows first, in ascending id order - the order the sender used
        received_ids = (numpy.concatenate(
            [numpy.intersect1d(gathered_owned[source], needed_ids)
             for source in range(self.size)]) if self.size
            else numpy.zeros(0, numpy.int64))
        if not numpy.array_equal(numpy.sort(received_ids), needed_ids):
            raise ValueError(
                "The halo exchange plan does not cover this rank's halo: some halo node "
                "is owned by no rank. The owned sets must partition the global nodes.")
        row_of = {int(node_id): row for row, node_id in enumerate(node_ids)}
        self.owned_rows = torch.tensor([row_of[int(i)] for i in owned_ids], dtype=torch.int64)
        self.halo_rows = torch.tensor([row_of[int(i)] for i in received_ids], dtype=torch.int64)

    def UnpadReceived(self, received):
        """The true rows of each source's padded block, concatenated."""
        torch = _TryImportTorch()
        if not self.block_size:
            return received
        keep = []
        for source in range(self.size):
            start = source * self.block_size
            keep.extend(range(start, start + self.true_sizes[source][self.rank]))
        return received[torch.tensor(keep, dtype=torch.int64, device=received.device)]


def ExchangeHaloFeatures(owned_features, plan: HaloExchangePlan, group=None):
    """The subgraph's full (N, C) feature tensor, differentiably.

    Args:
        owned_features: (n_owned, C) torch tensor, rows in ascending
            owned-node-id order - requires_grad is honoured.
        plan: A HaloExchangePlan for this rank.
        group: The torch process group (default: the world group).

    Returns:
        (N, C) tensor in the subgraph's node_ids order. Gradients reaching
        a halo row are summed back onto the owning rank's row by the
        collective's backward, so a loss over every rank's subgraph gives
        each owned row the total gradient of all its copies.
    """
    torch = _TryImportTorch()
    indexed_all_to_all_v = _TryImportAutogradCollectives()

    if owned_features.shape[0] != plan.n_owned:
        raise ValueError(
            f"owned_features has {owned_features.shape[0]} rows but this rank owns "
            f"{plan.n_owned} nodes; pass the owned rows in ascending id order.")
    full = owned_features.new_zeros((plan.n_nodes,) + tuple(owned_features.shape[1:]))
    full = full.index_copy(0, plan.owned_rows.to(full.device), owned_features)
    if not plan.halo_rows.numel() and not plan.block_size:
        return full

    padded = indexed_all_to_all_v(
        owned_features, plan.send_rows, plan.sizes, use_fp32=False, dim=0, group=group)
    received = plan.UnpadReceived(padded)
    if plan.halo_rows.numel():
        full = full.index_copy(0, plan.halo_rows.to(full.device), received)
    return full


def GatherOwnedPredictionsToRank0(model_part: Kratos.ModelPart, node_ids,
                                  owned_mask, values, data_communicator=None):
    """Assembles per-rank owned predictions into the serial layout on rank 0.

    Since owned sets partition the global node set exactly, concatenating
    them and sorting by node id reproduces what a single-rank run would have
    produced - which is what makes "distributed inference == serial
    inference" checkable.

    Returns:
        (ids, values) sorted by node id on rank 0; (None, None) elsewhere.
        Collective: call it on every rank.
    """
    if data_communicator is None:
        data_communicator = model_part.GetCommunicator().GetDataCommunicator()

    values = numpy.asarray(values, dtype=float)
    owned_mask = numpy.asarray(owned_mask, dtype=bool)
    owned_ids = numpy.asarray(node_ids, dtype=numpy.int64)[owned_mask]
    owned_values = values[owned_mask]
    width = int(owned_values.shape[1]) if owned_values.ndim > 1 else 1

    gathered_ids = data_communicator.GathervInts(
        [int(i) for i in owned_ids], 0)
    gathered_values = data_communicator.GathervDoubles(
        [float(v) for v in owned_values.reshape(-1)], 0)
    if data_communicator.Rank() != 0:
        return None, None

    ids = numpy.array([i for block in gathered_ids for i in block], dtype=numpy.int64)
    flat = numpy.array([v for block in gathered_values for v in block], dtype=float)
    stacked = flat.reshape(len(ids), width)
    order = numpy.argsort(ids)
    return ids[order], stacked[order]


def InitializeTorchProcessGroup(data_communicator=None, backend: str = "gloo",
                                address: str = "127.0.0.1", port: str = "29500"):
    """Initializes torch.distributed from Kratos's rank/size.

    Deliberately bypasses physicsnemo's DistributedManager: that helper
    passes a CUDA `device_id` to init_process_group, which requires an
    accelerator per rank and therefore fails on a single-GPU machine (and
    breaks gloo subgroup creation). Going direct to torch works with gloo on
    CPU, and coexists with a live Kratos MPI run in the same process.

    Returns:
        (rank, world_size). Idempotent: an already-initialized group is
        reused.
    """
    import os

    torch = _TryImportTorch()
    import torch.distributed as distributed

    if data_communicator is None:
        data_communicator = Kratos.ParallelEnvironment.GetDefaultDataCommunicator()
    rank, world_size = data_communicator.Rank(), data_communicator.Size()

    if not distributed.is_initialized():
        os.environ.setdefault("MASTER_ADDR", address)
        os.environ.setdefault("MASTER_PORT", port)
        os.environ["RANK"] = str(rank)
        os.environ["WORLD_SIZE"] = str(world_size)
        distributed.init_process_group(backend, rank=rank, world_size=world_size)
    return rank, world_size


def WrapForDataParallel(model):
    """Wraps a model in DistributedDataParallel for per-rank subgraphs.

    Note:
        DDP averages the per-rank losses, so with unequal node counts per
        rank the result is an unweighted-per-rank mean rather than a
        per-node one. Scale each rank's loss by its owned-node count if a
        per-node mean is what you want.

        Every rank must touch every parameter each step, or DDP needs
        find_unused_parameters=True.
    """
    torch = _TryImportTorch()
    from torch.nn.parallel import DistributedDataParallel

    return DistributedDataParallel(model)
