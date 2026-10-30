"""MPI tests for the matched process-group and device-mesh helpers.

physicsnemo's DistributedManager is initialized ONCE per MPI run (gloo, CPU)
with rank/size taken from the Kratos DataCommunicator: the launcher env
vars are popped so InitializeDistributedManager takes its deterministic
explicit-setup path instead of torch's auto-detection. Group names are
unique per test (physicsnemo forbids re-creating a name) and the manager is
never cleaned up between tests (that would destroy the shared process
group under the other ranks).
"""

import os

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.PhysicsNeMoApplication.distributed import distributed_utils
try:
    from physicsnemo.distributed.manager import DistributedManager
    have_physicsnemo = True
except ImportError:
    have_physicsnemo = False

_LAUNCHER_VARIABLES = ("RANK", "WORLD_SIZE", "LOCAL_RANK", "SLURM_PROCID", "OMPI_COMM_WORLD_RANK")


@KratosUnittest.skipUnless(have_physicsnemo, "Missing required python module: physicsnemo.")
class TestMpiProcessGroups(KratosUnittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._saved_environment = {
            name: os.environ.pop(name) for name in _LAUNCHER_VARIABLES if name in os.environ}
        distributed_utils.InitializeDistributedManager(backend="gloo", port="29601")

    @classmethod
    def tearDownClass(cls):
        os.environ.update(cls._saved_environment)

    def setUp(self):
        self.data_communicator = Kratos.Testing.GetDefaultDataCommunicator()

    def test_MatchedWorldGroup(self):
        world = self.data_communicator.Size()
        sub = distributed_utils.CreateMatchedProcessGroup(
            "pn_world", world, data_communicator=self.data_communicator)
        manager = DistributedManager()
        self.assertEqual(manager.group_size("pn_world"), world)
        self.assertEqual(manager.group_rank("pn_world"), self.data_communicator.Rank())
        self.assertEqual(sub.Size(), world)
        self.assertEqual(sub.Rank(), self.data_communicator.Rank())
        self.assertTrue(Kratos.ParallelEnvironment.HasDataCommunicator("pn_world"))

    def test_MatchedSizeOneGroups(self):
        sub = distributed_utils.CreateMatchedProcessGroup(
            "pn_singles", 1, data_communicator=self.data_communicator)
        manager = DistributedManager()
        self.assertEqual(manager.group_size("pn_singles"), 1)
        self.assertEqual(manager.group_rank("pn_singles"), 0)
        self.assertEqual(sub.Size(), 1)
        self.assertEqual(sub.Rank(), 0)

    def test_MatchedSizeTwoGroups(self):
        world = self.data_communicator.Size()
        if world % 2 != 0:
            self.skipTest(f"world size {world} is not divisible by 2")
        sub = distributed_utils.CreateMatchedProcessGroup(
            "pn_pairs", 2, data_communicator=self.data_communicator)
        self.assertEqual(sub.Size(), 2)
        self.assertEqual(sub.Rank(), self.data_communicator.Rank() % 2)

    def test_ParametersDrivenGroups(self):
        groups = distributed_utils.CreateMatchedProcessGroups(Kratos.Parameters("""{
            "process_groups" : [ { "name" : "pn_parameters", "size" : 1 } ]
        }"""), data_communicator=self.data_communicator)
        self.assertEqual(set(groups), {"pn_parameters"})
        self.assertEqual(groups["pn_parameters"].Size(), 1)

    def test_DeviceMeshValidationAndCreation(self):
        world = self.data_communicator.Size()
        with self.assertRaisesRegex(ValueError, "implies"):
            distributed_utils.InitializeDeviceMesh(
                (world + 1,), ("data",), data_communicator=self.data_communicator)
        with self.assertRaisesRegex(ValueError, "entries but"):
            distributed_utils.InitializeDeviceMesh(
                (world, 1), ("data",), data_communicator=self.data_communicator)
        mesh = distributed_utils.InitializeDeviceMesh(
            (-1,), ("data",), data_communicator=self.data_communicator)
        self.assertIsNotNone(mesh)


if __name__ == '__main__':
    KratosUnittest.main()
