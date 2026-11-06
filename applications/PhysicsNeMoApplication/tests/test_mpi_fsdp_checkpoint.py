"""Sharded FSDP2 checkpoints across ranks.

Run with:
    mpiexec -np 2 python3 test_PhysicsNeMoApplication_mpi.py

Everything here is gloo on CPU. That proves the sharding and gathering
logic but NOT the accelerator transport: NCCL rejects two ranks on one
GPU, so a multi-GPU NCCL round trip remains untested on this machine -
the same honest limit TestMpiDataParallelTraining ships under.
"""

import io
import os
import zipfile
from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

from KratosMultiphysics.PhysicsNeMoApplication.training import training_utils
try:
    import torch
    import torch.distributed as distributed
    from torch.distributed.fsdp import fully_shard
    from torch.distributed.tensor import DTensor
    have_fsdp = True
except ImportError:
    have_fsdp = False

try:
    import physicsnemo  # noqa: F401
    from physicsnemo.models.mlp.fully_connected import FullyConnected  # noqa: F401
    have_physicsnemo = True
except ImportError:
    have_physicsnemo = False


@KratosUnittest.skipUnless(have_fsdp, "Missing torch with FSDP2.")
@KratosUnittest.skipUnless(have_physicsnemo, "Missing required python module: physicsnemo.")
class TestMpiFsdpCheckpoint(KratosUnittest.TestCase):
    """A genuinely sharded save/load round trip."""

    _PORT = "29602"

    @classmethod
    def setUpClass(cls):
        from torch.distributed.device_mesh import init_device_mesh
        from KratosMultiphysics.PhysicsNeMoApplication.distributed import graph_partition_utils
        cls.data_communicator = Kratos.Testing.GetDefaultDataCommunicator()
        if not distributed.is_initialized():
            # InitializeTorchProcessGroup uses setdefault for MASTER_PORT, so
            # a port left in the environment by an earlier test class in the
            # same MPI run would win - and rendezvous on it fails once that
            # test's group is gone. Set it explicitly.
            os.environ["MASTER_ADDR"] = "127.0.0.1"
            os.environ["MASTER_PORT"] = cls._PORT
            graph_partition_utils.InitializeTorchProcessGroup(
                cls.data_communicator, backend="gloo", port=cls._PORT)
        # Built ONCE. Creating a device mesh per test builds fresh
        # sub-groups each time, which deadlocks when another test class has
        # already brought the default group up (test_mpi_distributed_groups
        # initializes one and deliberately never tears it down).
        cls.mesh = init_device_mesh("cpu", (distributed.get_world_size(),),
                                    mesh_dim_names=("dp",))

    @classmethod
    def tearDownClass(cls):
        # Deliberately NOT destroying the process group: it may be shared
        # with other test classes in the same MPI run, and tearing it down
        # here would strand them. Same convention as
        # test_mpi_distributed_groups.
        pass

    def setUp(self):
        self.checkpoint = Path("test_mpi_fsdp_checkpoint.mdlus")

    def tearDown(self):
        self.data_communicator.Barrier()
        if self.data_communicator.Rank() == 0:
            KratosUtilities.DeleteFileIfExisting(str(self.checkpoint))

    @staticmethod
    def _Model():
        from physicsnemo.models.mlp.fully_connected import FullyConnected
        torch.manual_seed(0)     # identical initial weights on every rank
        return FullyConnected(in_features=4, out_features=2,
                              layer_size=16, num_layers=2)

    @staticmethod
    def _Shard(model):
        """fully_shard onto an EXPLICIT CPU mesh.

        Without a mesh argument, fully_shard builds a CUDA mesh whenever a
        CUDA device is visible - even though the process group here is gloo
        on CPU. Gathering those DTensors then segfaults inside gloo's
        allgather (measured on torch 2.13). Naming the CPU mesh is both the
        fix and the honest thing to do for a CPU/gloo run.
        """
        fully_shard(model, mesh=TestMpiFsdpCheckpoint.mesh)
        return model

    def test_ParametersAreTrulySharded(self):
        # the precondition everything else rests on: if local == global the
        # "sharded" round trip below would prove nothing
        model = self._Shard(self._Model())
        sharded = 0
        for parameter in model.parameters():
            if isinstance(parameter, DTensor):
                if tuple(parameter.to_local().shape) != tuple(parameter.shape):
                    sharded += 1
        self.assertGreater(sharded, 0,
                           msg="no parameter was actually split across ranks")

    def test_ShardedSaveGathersAndReloads(self):
        model = self._Model()
        reference = {name: value.detach().clone()
                     for name, value in model.state_dict().items()}
        self._Shard(model)
        self.assertTrue(any(isinstance(p, DTensor) for p in model.parameters()))

        # collective: every rank must call it, only rank 0 writes
        checkpoint_type = training_utils.SaveTrainedModel(model, str(self.checkpoint))
        self.assertEqual(checkpoint_type, "physicsnemo")
        self.data_communicator.Barrier()

        if self.data_communicator.Rank() == 0:
            self.assertTrue(self.checkpoint.is_file())
            with zipfile.ZipFile(str(self.checkpoint)) as archive:
                stored = torch.load(io.BytesIO(archive.read("model.pt")),
                                    weights_only=False)
            # the bug this guards: DTensors written straight to disk
            self.assertEqual({type(v).__name__ for v in stored.values()}, {"Tensor"})
            for name, value in stored.items():
                self.assertEqual(tuple(value.shape), tuple(reference[name].shape),
                                 msg=f"{name} was written as a local shard")

            reloaded = physicsnemo.Module.from_checkpoint(str(self.checkpoint))
            for name, value in reloaded.state_dict().items():
                self.assertTrue(torch.allclose(reference[name], value),
                                msg=f"{name} differs from the pre-shard weights")
        else:
            # the ranks share a working directory, so this rank can see the
            # file rank 0 wrote; what matters is that it holds the GLOBAL
            # tensors rather than having been clobbered with a local shard
            with zipfile.ZipFile(str(self.checkpoint)) as archive:
                stored = torch.load(io.BytesIO(archive.read("model.pt")),
                                    weights_only=False)
            for name, value in stored.items():
                self.assertEqual(tuple(value.shape), tuple(reference[name].shape))

    def test_OnlyRankZeroWrites(self):
        # with a rank-distinct path the "one writer" contract is unambiguous:
        # before the fix every rank wrote its own local shard
        model = self._Shard(self._Model())
        rank = self.data_communicator.Rank()
        path = Path(f"test_mpi_fsdp_rank_{rank}.mdlus")
        try:
            training_utils.SaveTrainedModel(model, str(path))
            self.data_communicator.Barrier()
            self.assertEqual(path.is_file(), rank == 0)
        finally:
            KratosUtilities.DeleteFileIfExisting(str(path))


if __name__ == '__main__':
    KratosUnittest.main()
