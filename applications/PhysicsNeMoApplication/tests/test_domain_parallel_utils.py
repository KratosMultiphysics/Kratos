"""Serial (world_size = 1) coverage of domain_parallel_utils.

The sharding itself is asserted across ranks in the MPI suite
(test_mpi_domain_parallel); this pins the single-rank contract - every
helper is a valid no-op-shaped call at world size 1 - and the refusals,
over a one-process gloo group (the pattern of
test_training_utils.TestSaveShardedModel).
"""

import os
from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest

try:
    import torch
    import torch.distributed as distributed
    have_torch = True
except ImportError:
    have_torch = False

try:
    from physicsnemo.domain_parallel import ShardTensor  # noqa: F401
    have_shard_tensor = have_torch
except ImportError:
    have_shard_tensor = False


@KratosUnittest.skipUnless(have_shard_tensor,
                           "Missing required python modules: torch (>= 2.6), physicsnemo.domain_parallel.")
class TestDomainParallelUtilsSerial(KratosUnittest.TestCase):
    _PORT = "29672"

    def setUp(self):
        from KratosMultiphysics.PhysicsNeMoApplication.distributed import domain_parallel_utils
        self.utils = domain_parallel_utils
        self.previous = {key: os.environ.get(key)
                         for key in ("MASTER_ADDR", "MASTER_PORT", "RANK", "WORLD_SIZE")}
        os.environ.update(MASTER_ADDR="127.0.0.1", MASTER_PORT=self._PORT,
                          RANK="0", WORLD_SIZE="1")
        self.started_group = not distributed.is_initialized()
        if self.started_group:
            distributed.init_process_group("gloo", rank=0, world_size=1)
        self.mesh = self.utils.InitializeDomainMesh(port=self._PORT)

        self.model = Kratos.Model()
        self.model_part = self.model.CreateModelPart("Main")
        self.model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        self.model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        for i in range(5):
            node = self.model_part.CreateNewNode(10 * (i + 1), float(i), 2.0 * i, 0.0)
            node.SetSolutionStepValue(Kratos.PRESSURE, 1.0 + i)
            node.SetSolutionStepValue(Kratos.VELOCITY, [float(i), -float(i), 0.5 * i])

    def tearDown(self):
        if self.started_group and distributed.is_initialized():
            distributed.destroy_process_group()
        for key, value in self.previous.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value

    def test_MeshIsOneCpuRank(self):
        self.assertEqual(self.mesh.size(), 1)
        self.assertEqual(self.mesh.device_type, "cpu")
        with self.assertRaisesRegex(ValueError, "device_type"):
            self.utils.InitializeDomainMesh(device_type="tpu")

    def test_ShardLocalRowsRoundTrips(self):
        from torch.distributed.tensor import Shard
        rows = numpy.arange(15.0).reshape(5, 3)
        shard = self.utils.ShardLocalRows(rows, self.mesh)
        self.assertEqual(tuple(shard.placements), (Shard(0),))
        numpy.testing.assert_array_equal(shard.full_tensor().numpy(), rows)
        # a negative axis is normalized
        along_last = self.utils.ShardLocalRows(torch.from_numpy(rows), self.mesh, dim=-1)
        self.assertEqual(tuple(along_last.placements), (Shard(1),))
        with self.assertRaisesRegex(ValueError, "at least one axis"):
            self.utils.ShardLocalRows(numpy.float64(1.0), self.mesh)

    def test_ShardKratosFieldThenWriteBack(self):
        specs = [("PRESSURE", "node_historical"), ("VELOCITY", "node_historical")]
        shard, owned_ids = self.utils.ShardKratosField(self.model_part, specs, self.mesh)
        numpy.testing.assert_array_equal(owned_ids, [10, 20, 30, 40, 50])
        expected = numpy.array([[1.0 + i, float(i), -float(i), 0.5 * i] for i in range(5)])
        numpy.testing.assert_allclose(shard.full_tensor().numpy(), expected, rtol=1e-14)
        with self.assertRaisesRegex(ValueError, "at least one field"):
            self.utils.ShardKratosField(self.model_part, [], self.mesh)

        self.utils.WriteShardedField(shard * 2.0 + 1.0, self.model_part,
                                     [("TEMPERATURE", "node_non_historical"),
                                      ("DISPLACEMENT", "node_non_historical")])
        for i, node in enumerate(self.model_part.Nodes):
            self.assertAlmostEqual(node.GetValue(Kratos.TEMPERATURE), 2.0 * (1.0 + i) + 1.0)
            numpy.testing.assert_allclose(node.GetValue(Kratos.DISPLACEMENT),
                                          [2.0 * i + 1.0, -2.0 * i + 1.0, i + 1.0])
        # the card normalization reaches this write path too
        self.utils.WriteShardedField(
            shard, self.model_part, [("TEMPERATURE", "node_non_historical"),
                                     ("DISPLACEMENT", "node_non_historical")],
            normalization={"type": "mean_std", "mean": [100.0], "std": [2.0]})
        self.assertAlmostEqual(self.model_part.GetNode(10).GetValue(Kratos.TEMPERATURE), 102.0)

    def test_ShardGridAlongAxis(self):
        grid = numpy.random.default_rng(0).standard_normal((1, 2, 4, 6))
        sharded = self.utils.ShardGridAlongAxis(grid, self.mesh, axis=3)
        numpy.testing.assert_array_equal(sharded.full_tensor().numpy(), grid)
        with self.assertRaisesRegex(ValueError, "channel axis"):
            self.utils.ShardGridAlongAxis(grid, self.mesh, axis=0)
        with self.assertRaisesRegex(ValueError, "fewer than"):
            self.utils.ShardGridAlongAxis(numpy.zeros((2, 0)), self.mesh, axis=1)

    def test_MeshWideValueOfAReduction(self):
        rows = numpy.arange(12.0).reshape(4, 3)
        shard = self.utils.ShardLocalRows(rows, self.mesh)
        self.assertAlmostEqual(self.utils.MeshWideValue(shard.square().mean()),
                               float(numpy.square(rows).mean()), places=12)
        self.assertAlmostEqual(self.utils.MeshWideValue(torch.tensor(2.5)), 2.5)
        self.assertAlmostEqual(self.utils.MeshWideValue(2.5), 2.5)


if __name__ == '__main__':
    KratosUnittest.main()
