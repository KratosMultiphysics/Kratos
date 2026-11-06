"""Domain parallelism over ShardTensor across MPI ranks.

Run with:
    mpiexec -np 2 python3 test_PhysicsNeMoApplication_mpi.py

Everything here is gloo on CPU, with no DistributedManager: that proves the
sharding, the halo exchange and the gradient mathematics but NOT the
accelerator transport - NCCL rejects two ranks on one GPU - the same honest
limit TestMpiFsdpCheckpoint and TestMpiDataParallelTraining ship under. The
roadmap used to call this path hardware-blocked outright; the actual
CUDA-only remnants are listed in domain_parallel_utils' docstring.
"""

import copy
import os

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest

from test_mpi_graph_partition import _CreateDistributedModelPart, _CreateSerialReference

try:
    import torch
    import torch.distributed as distributed
    have_torch = True
except ImportError:
    have_torch = False

try:
    from physicsnemo.domain_parallel import ShardTensor  # noqa: F401
    from physicsnemo.models.mlp.fully_connected import FullyConnected  # noqa: F401
    have_shard_tensor = have_torch
except ImportError:
    have_shard_tensor = False

_FIELD = [("PRESSURE", "node_historical")]


@KratosUnittest.skipUnless(have_shard_tensor,
                           "Missing required python modules: torch (>= 2.6), physicsnemo.domain_parallel.")
class TestMpiDomainParallel(KratosUnittest.TestCase):
    """One Kratos field sharded across the ranks, exercised through the
    ops physicsnemo makes mesh-aware, always against the serial answer."""

    _PORT = "29603"

    @classmethod
    def setUpClass(cls):
        from KratosMultiphysics.PhysicsNeMoApplication.distributed import domain_parallel_utils
        cls.utils = domain_parallel_utils
        cls.data_communicator = Kratos.Testing.GetDefaultDataCommunicator()
        if not distributed.is_initialized():
            # InitializeTorchProcessGroup uses setdefault for MASTER_PORT; a
            # port left by an earlier test class would win (see
            # test_mpi_fsdp_checkpoint). Set it explicitly.
            os.environ["MASTER_ADDR"] = "127.0.0.1"
            os.environ["MASTER_PORT"] = cls._PORT
        # built ONCE: a mesh creates sub-groups, and doing that per test
        # while other classes' groups are live deadlocks
        cls.mesh = cls.utils.InitializeDomainMesh(cls.data_communicator, port=cls._PORT)

    @classmethod
    def tearDownClass(cls):
        pass  # the process group is shared with the other MPI test classes

    def setUp(self):
        self.model = Kratos.Model()
        self.model_part, _ = _CreateDistributedModelPart(self.model)
        self.reference = _CreateSerialReference(self.model)
        self.world_size = self.data_communicator.Size()

    def _SerialField(self):
        ids = numpy.array([node.Id for node in self.reference.Nodes], dtype=numpy.int64)
        values = numpy.array([node.GetSolutionStepValue(Kratos.PRESSURE)
                              for node in self.reference.Nodes])
        order = numpy.argsort(ids)
        return ids[order], values[order]

    def _GatherIds(self, owned_ids):
        pieces = [None] * self.world_size
        distributed.all_gather_object(pieces, owned_ids.tolist())
        return numpy.concatenate([numpy.asarray(p, dtype=numpy.int64) for p in pieces])

    def test_ShardedFieldGathersToTheSerialLayout(self):
        shard, owned_ids = self.utils.ShardKratosField(self.model_part, _FIELD, self.mesh)
        # the precondition everything else rests on: the ranks really hold
        # different, uneven-or-not, PIECES of one tensor
        self.assertLess(shard.to_local().shape[0], shard.shape[0])
        self.assertEqual(shard.to_local().shape[0], len(owned_ids))

        ids = self._GatherIds(owned_ids)
        values = shard.full_tensor().numpy()[:, 0]
        self.assertEqual(len(ids), self.reference.NumberOfNodes())
        order = numpy.argsort(ids)
        serial_ids, serial_values = self._SerialField()
        numpy.testing.assert_array_equal(ids[order], serial_ids)
        numpy.testing.assert_allclose(values[order], serial_values, rtol=1e-14)

    def test_PointwiseModelOnTheShardEqualsSerial(self):
        from physicsnemo.models.mlp.fully_connected import FullyConnected
        shard, _ = self.utils.ShardKratosField(self.model_part, _FIELD, self.mesh)
        torch.manual_seed(0)
        model = FullyConnected(in_features=1, out_features=2, layer_size=8, num_layers=2).double()

        output = model(shard)
        self.assertEqual(output.to_local().shape[0], shard.to_local().shape[0])
        with torch.no_grad():
            serial = model(shard.full_tensor())
        numpy.testing.assert_allclose(output.full_tensor().detach().numpy(), serial.numpy(),
                                      rtol=1e-12)

    def test_HaloConvolutionOnAShardedGridEqualsSerial(self):
        # a (1, C, H, W) grid split along W: each rank's 3x3 convolution needs
        # one column from each neighbour, which is the halo exchange
        y, x = numpy.meshgrid(numpy.linspace(0.0, 1.0, 8), numpy.linspace(0.0, 1.0, 16),
                              indexing="ij")
        grid = numpy.stack([x + 2.0 * y, numpy.sin(3.0 * x) * numpy.cos(2.0 * y)])[None]
        sharded = self.utils.ShardGridAlongAxis(grid, self.mesh, axis=3)
        self.assertLess(sharded.to_local().shape[3], grid.shape[3])

        torch.manual_seed(1)
        conv = torch.nn.Conv2d(2, 3, kernel_size=3, padding=1).double()
        with torch.no_grad():
            serial = conv(torch.from_numpy(grid))
        result = conv(sharded).full_tensor().detach()
        numpy.testing.assert_allclose(result.numpy(), serial.numpy(), rtol=1e-12, atol=1e-14)

    def test_GradientOfAMeshWideLossIsTheSerialGradientOnEveryRank(self):
        """The training step. A mean over a ShardTensor is a mesh-wide
        quantity whose backward already delivers the FULL gradient on every
        rank - so an extra DDP-style all-reduce would be wrong by the rank
        count. Both halves are pinned."""
        from physicsnemo.models.mlp.fully_connected import FullyConnected
        shard, _ = self.utils.ShardKratosField(self.model_part, _FIELD, self.mesh)
        torch.manual_seed(2)
        model = FullyConnected(in_features=1, out_features=2, layer_size=8, num_layers=2).double()
        serial_model = copy.deepcopy(model)

        loss = (model(shard) - 1.0).square().mean()
        serial_loss = (serial_model(shard.full_tensor()) - 1.0).square().mean()
        self.assertAlmostEqual(self.utils.MeshWideValue(loss), float(serial_loss.detach()),
                               places=12)
        loss.backward()
        serial_loss.backward()

        for (name, parameter), serial_parameter in zip(model.named_parameters(),
                                                       serial_model.parameters()):
            with self.subTest(parameter=name):
                gradient = parameter.grad.detach().clone()
                numpy.testing.assert_allclose(gradient.numpy(), serial_parameter.grad.numpy(),
                                              rtol=1e-6)
                # identical on every rank
                reference = gradient.clone()
                distributed.broadcast(reference, src=0)
                numpy.testing.assert_allclose(gradient.numpy(), reference.numpy(), rtol=1e-9)
                # and summing again over the ranks would NOT be the gradient
                summed = gradient.clone()
                distributed.all_reduce(summed, op=distributed.ReduceOp.SUM)
                self.assertFalse(numpy.allclose(summed.numpy(), serial_parameter.grad.numpy(),
                                                rtol=1e-3))

    def test_WriteBackReachesOwnedRowsAndGhosts(self):
        shard, _ = self.utils.ShardKratosField(self.model_part, _FIELD, self.mesh)
        # a GLOBAL precondition: in the slab fixture the interface column is
        # owned by the left rank, so rank 0 holds no ghosts at all - and a
        # rank-local assertion here would skip the collective ghost sync on
        # one rank and hang the others (it did)
        ghost_nodes = self.data_communicator.SumAll(
            int(self.model_part.GetCommunicator().GhostMesh().NumberOfNodes()))
        self.assertGreater(ghost_nodes, 0)
        self.utils.WriteShardedField(shard * 2.0 + 1.0, self.model_part,
                                     [("TEMPERATURE", "node_non_historical")])
        for node in self.model_part.Nodes:  # owned AND ghost
            self.assertAlmostEqual(node.GetValue(Kratos.TEMPERATURE),
                                   2.0 * (node.X + 2.0 * node.Y) + 1.0, places=12)
        # with a card normalization on the write path, like every other writer
        self.utils.WriteShardedField(
            shard, self.model_part, [("TEMPERATURE", "node_non_historical")],
            normalization={"type": "mean_std", "mean": [10.0], "std": [3.0]})
        for node in self.model_part.Nodes:
            self.assertAlmostEqual(node.GetValue(Kratos.TEMPERATURE),
                                   3.0 * (node.X + 2.0 * node.Y) + 10.0, places=12)


if __name__ == '__main__':
    KratosUnittest.main()
