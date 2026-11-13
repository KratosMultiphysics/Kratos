"""Distributed graph partitioning: per-rank halo subgraphs and DDP training.

Run with:
    mpiexec -np 2 python3 test_PhysicsNeMoApplication_mpi.py
"""

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.PhysicsNeMoApplication.bridges import graph_bridge
from KratosMultiphysics.PhysicsNeMoApplication.distributed import graph_partition_utils
try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

try:
    from physicsnemo.nn.module.gnn_layers.graph_types import PYG_AVAILABLE
    from physicsnemo.models.meshgraphnet import MeshGraphNet
    have_meshgraphnet = have_torch and PYG_AVAILABLE
except ImportError:
    have_meshgraphnet = False


def _CreateDistributedModelPart(model, divisions=6):
    """The Metis-free slab fixture: each rank builds its own x-slab."""
    import KratosMultiphysics.mpi as KratosMPI

    data_communicator = Kratos.Testing.GetDefaultDataCommunicator()
    rank, size = data_communicator.Rank(), data_communicator.Size()

    model_part = model.CreateModelPart("Partitioned")
    model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 2
    model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
    model_part.AddNodalSolutionStepVariable(Kratos.PARTITION_INDEX)
    properties = model_part.CreateNewProperties(1)
    n = divisions

    def cell_owner(cell_x):
        return min(cell_x * size // n, size - 1)

    def node_owner(node_x):
        return cell_owner(max(node_x - 1, 0))

    def node_id(node_x, node_y):
        return node_x * (n + 1) + node_y + 1

    my_columns = [cell_x for cell_x in range(n) if cell_owner(cell_x) == rank]
    needed = set()
    for cell_x in my_columns:
        for cell_y in range(n):
            for dx, dy in ((0, 0), (1, 0), (1, 1), (0, 1)):
                needed.add((cell_x + dx, cell_y + dy))
    for node_x, node_y in sorted(needed):
        node = model_part.CreateNewNode(node_id(node_x, node_y), node_x / n, node_y / n, 0.0)
        node.SetSolutionStepValue(Kratos.PARTITION_INDEX, node_owner(node_x))
    for cell_x in my_columns:
        for cell_y in range(n):
            quad = [node_id(cell_x, cell_y), node_id(cell_x + 1, cell_y),
                    node_id(cell_x + 1, cell_y + 1), node_id(cell_x, cell_y + 1)]
            base = 2 * (cell_x * n + cell_y)
            model_part.CreateNewElement("Element2D3N", base + 1, quad[:3], properties)
            model_part.CreateNewElement("Element2D3N", base + 2, [quad[0], quad[2], quad[3]], properties)

    KratosMPI.ParallelFillCommunicator(model_part, data_communicator).Execute()
    for node in model_part.Nodes:
        node.SetSolutionStepValue(Kratos.PRESSURE, node.X + 2.0 * node.Y)
    return model_part, data_communicator


def _CreateSerialReference(model, divisions=6):
    """The same mesh on one rank, as the ground truth."""
    model_part = model.CreateModelPart("SerialReference")
    model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 2
    model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
    properties = model_part.CreateNewProperties(1)
    n = divisions

    def node_id(node_x, node_y):
        return node_x * (n + 1) + node_y + 1

    for node_x in range(n + 1):
        for node_y in range(n + 1):
            model_part.CreateNewNode(node_id(node_x, node_y), node_x / n, node_y / n, 0.0)
    for cell_x in range(n):
        for cell_y in range(n):
            quad = [node_id(cell_x, cell_y), node_id(cell_x + 1, cell_y),
                    node_id(cell_x + 1, cell_y + 1), node_id(cell_x, cell_y + 1)]
            base = 2 * (cell_x * n + cell_y)
            model_part.CreateNewElement("Element2D3N", base + 1, quad[:3], properties)
            model_part.CreateNewElement("Element2D3N", base + 2, [quad[0], quad[2], quad[3]], properties)
    for node in model_part.Nodes:
        node.SetSolutionStepValue(Kratos.PRESSURE, node.X + 2.0 * node.Y)
    return model_part


def _OneHopSum(edge_index, values):
    """y_i = sum over neighbours j of x_j - the message-passing primitive."""
    result = numpy.zeros_like(values)
    for sender, receiver in zip(edge_index[0], edge_index[1]):
        result[receiver] += values[sender]
    return result


class TestMpiHaloSubgraph(KratosUnittest.TestCase):

    def setUp(self):
        self.model = Kratos.Model()
        self.model_part, self.data_communicator = _CreateDistributedModelPart(self.model)

    def test_OwnedSetsPartitionTheGlobalNodeSetExactly(self):
        _, _, _, node_ids, owned_mask = graph_partition_utils.BuildHaloSubgraph(
            self.model_part)
        owned_count = int(owned_mask.sum())
        total = self.data_communicator.SumAll(owned_count)
        self.assertEqual(total, self.model_part.GetCommunicator().GlobalNumberOfNodes())
        # halo rows exist, and they are not this rank's to predict
        self.assertGreaterEqual(len(node_ids), owned_count)

    def test_OneHopMatchesSerialAtEveryOwnedNode(self):
        """The property the halo exists for."""
        reference_model = Kratos.Model()
        reference_part = _CreateSerialReference(reference_model)
        _, reference_edges, _, reference_ids = graph_bridge.BuildGraph(reference_part)
        reference_values = numpy.array(
            [node.GetSolutionStepValue(Kratos.PRESSURE) for node in reference_part.Nodes])
        reference_sum = _OneHopSum(reference_edges, reference_values)
        reference_of = {int(i): r for r, i in enumerate(reference_ids)}

        features, edge_index, _, node_ids, owned_mask = \
            graph_partition_utils.BuildHaloSubgraph(
                self.model_part, field_specs=(("PRESSURE", "node_historical"),))
        local_sum = _OneHopSum(edge_index, features[:, 0])

        mismatches = [int(node_ids[row]) for row in range(len(node_ids))
                      if owned_mask[row]
                      and abs(local_sum[row] - reference_sum[reference_of[int(node_ids[row])]]) > 1e-9]
        self.assertEqual(mismatches, [], msg=f"halo subgraph wrong at owned nodes {mismatches}")

    def test_PlainBuildGraphIsWrongAtTheInterface(self):
        """Negative control: the bug the halo subgraph fixes.

        BuildGraph on a distributed part truncates neighbourhoods at the
        partition boundary, and does so even at nodes this rank OWNS.
        """
        reference_model = Kratos.Model()
        reference_part = _CreateSerialReference(reference_model)
        _, reference_edges, _, reference_ids = graph_bridge.BuildGraph(reference_part)
        reference_values = numpy.array(
            [node.GetSolutionStepValue(Kratos.PRESSURE) for node in reference_part.Nodes])
        reference_sum = _OneHopSum(reference_edges, reference_values)
        reference_of = {int(i): r for r, i in enumerate(reference_ids)}

        features, edge_index, _, node_ids = graph_bridge.BuildGraph(
            self.model_part, (("PRESSURE", "node_historical"),))
        local_sum = _OneHopSum(edge_index, features[:, 0])
        mismatches = [int(node_ids[row]) for row in range(len(node_ids))
                      if abs(local_sum[row] - reference_sum[reference_of[int(node_ids[row])]]) > 1e-9]
        # with more than one rank there is always a truncated interface
        if self.data_communicator.Size() > 1:
            self.assertGreater(len(mismatches), 0)

    def test_MoreRingsGrowTheHalo(self):
        _, _, _, ids_one, owned_one = graph_partition_utils.BuildHaloSubgraph(
            self.model_part, num_halo_rings=1)
        _, _, _, ids_two, owned_two = graph_partition_utils.BuildHaloSubgraph(
            self.model_part, num_halo_rings=2)
        if self.data_communicator.Size() > 1:
            self.assertGreaterEqual(len(ids_two), len(ids_one))
        # the owned set never changes with halo depth
        self.assertEqual(int(owned_one.sum()), int(owned_two.sum()))

    def test_GatheredOwnedValuesReproduceTheSerialLayout(self):
        features, _, _, node_ids, owned_mask = graph_partition_utils.BuildHaloSubgraph(
            self.model_part, field_specs=(("PRESSURE", "node_historical"),))
        ids, values = graph_partition_utils.GatherOwnedPredictionsToRank0(
            self.model_part, node_ids, owned_mask, features)
        if self.data_communicator.Rank() != 0:
            return
        reference_model = Kratos.Model()
        reference_part = _CreateSerialReference(reference_model)
        expected = numpy.array(
            [node.GetSolutionStepValue(Kratos.PRESSURE) for node in reference_part.Nodes])
        self.assertEqual(len(ids), len(expected))
        numpy.testing.assert_allclose(values[:, 0], expected, atol=1e-12)


@KratosUnittest.skipUnless(have_meshgraphnet,
                           "Missing required python modules: physicsnemo with torch_geometric/torch_scatter.")
class TestMpiDataParallelTraining(KratosUnittest.TestCase):

    def setUp(self):
        self.model = Kratos.Model()
        self.model_part, self.data_communicator = _CreateDistributedModelPart(self.model)

    def test_GradientsAreIdenticalAcrossRanksAfterDdpBackward(self):
        """Data-parallel GNN training over per-rank subgraphs.

        Uses gloo directly rather than physicsnemo's DistributedManager,
        which requires a CUDA device per rank and so cannot run here.
        """
        import torch.distributed as distributed

        rank, world = graph_partition_utils.InitializeTorchProcessGroup(
            self.data_communicator, port="29601")
        try:
            features, edge_index, edge_features, node_ids, owned_mask = \
                graph_partition_utils.BuildHaloSubgraph(
                    self.model_part, field_specs=(("PRESSURE", "node_historical"),))

            torch.manual_seed(0)   # identical initial weights on every rank
            model = MeshGraphNet(input_dim_nodes=1, input_dim_edges=4, output_dim=1,
                                 processor_size=2, hidden_dim_processor=8,
                                 hidden_dim_node_encoder=8, hidden_dim_edge_encoder=8,
                                 hidden_dim_node_decoder=8).double()
            wrapped = graph_partition_utils.WrapForDataParallel(model)

            graph = graph_bridge.ToPyGGraph(edge_index, len(node_ids))
            prediction = wrapped(torch.from_numpy(features),
                                 torch.from_numpy(edge_features), graph)
            # score only what this rank owns
            owned = torch.from_numpy(owned_mask)
            loss = prediction[owned].square().mean()
            loss.backward()

            gradients = [p.grad for p in model.parameters() if p.grad is not None]
            self.assertGreater(len(gradients), 0)
            for gradient in gradients:
                reference = gradient.clone()
                distributed.broadcast(reference, src=0)
                self.assertTrue(torch.equal(gradient, reference),
                                msg="DDP did not allreduce the gradients")
        finally:
            if distributed.is_initialized():
                distributed.destroy_process_group()


@KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
class TestMpiDifferentiableHaloExchange(KratosUnittest.TestCase):
    """The halo exchange through torch autograd: the same values the eager
    exchange copies in, plus a gradient path back to the owning rank."""

    def setUp(self):
        self.model = Kratos.Model()
        self.model_part, self.data_communicator = _CreateDistributedModelPart(self.model)
        graph_partition_utils.InitializeTorchProcessGroup(self.data_communicator)
        (self.node_features, _, _, self.node_ids,
         self.owned_mask) = graph_partition_utils.BuildHaloSubgraph(
            self.model_part, 1, field_specs=[("PRESSURE", "node_historical")])
        self.plan = graph_partition_utils.HaloExchangePlan(
            self.node_ids, self.owned_mask, self.data_communicator)

    def test_ExchangedFeaturesMatchTheEagerHalo(self):
        owned = torch.tensor(self.node_features[self.owned_mask], dtype=torch.float64)
        full = graph_partition_utils.ExchangeHaloFeatures(owned, self.plan)
        numpy.testing.assert_array_equal(full.detach().numpy(), self.node_features)

    def test_TheGradientOfEveryCopyReturnsToTheOwner(self):
        """Weight rank r's whole subgraph by r + 1 and sum over ranks. An
        owned row's gradient is then its own weight plus the weight of every
        rank holding a halo copy of it - the sum the collective's backward
        must deliver to the owner."""
        rank = self.data_communicator.Rank()
        owned = torch.tensor(
            self.node_features[self.owned_mask], dtype=torch.float64, requires_grad=True)
        full = graph_partition_utils.ExchangeHaloFeatures(owned, self.plan)
        (full * float(rank + 1)).sum().backward()

        halo_ids = [int(i) for i in self.node_ids[~self.owned_mask]]
        gathered_halos = self.data_communicator.AllGathervInts(halo_ids)
        owned_ids = numpy.sort(self.node_ids[self.owned_mask])
        expected = numpy.full(len(owned_ids), float(rank + 1))
        for other, ids in enumerate(gathered_halos):
            copies = set(int(i) for i in ids)
            for row, node_id in enumerate(owned_ids):
                if int(node_id) in copies:
                    expected[row] += float(other + 1)
        numpy.testing.assert_allclose(owned.grad.numpy().ravel(), expected, rtol=0, atol=1e-12)

        # somewhere an owned node really is someone else's halo - summed
        # over ALL ranks, so no rank skips the collective
        interface_rows = int((expected > rank + 1).sum())
        if self.data_communicator.Size() > 1:
            self.assertGreater(self.data_communicator.SumAll(interface_rows), 0)

    def test_TheExchangeIsRaggedAndIsPaddedForGloo(self):
        """The reason the plan pads at all: a halo exchange is ragged by
        nature - a rank sends rows to the neighbour it shares an interface
        with and none to itself - and gloo's all-to-all rejects blocks of
        differing shape. Every block is padded to one size and sliced back."""
        rank = self.data_communicator.Rank()
        if self.data_communicator.Size() > 1:
            counts = self.plan.true_sizes[rank]
            self.assertGreater(max(counts), 0)
            self.assertEqual(counts[rank], 0)      # a rank needs nothing of its own
            self.assertNotEqual(min(counts), max(counts))   # genuinely ragged
        self.assertTrue(all(size == self.plan.block_size
                            for row in self.plan.sizes for size in row))
        self.assertTrue(all(rows.numel() == self.plan.block_size
                            for rows in self.plan.send_rows))

    def test_AWrongOwnedRowCountIsRefused(self):
        """Refused before any collective, identically on every rank, so the
        error cannot leave a peer waiting."""
        wrong = torch.zeros(self.plan.n_owned + 1, 1, dtype=torch.float64)
        with self.assertRaisesRegex(ValueError, "owned"):
            graph_partition_utils.ExchangeHaloFeatures(wrong, self.plan)


if __name__ == '__main__':
    KratosUnittest.main()
