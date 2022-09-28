import KratosMultiphysics as Kratos

import KratosMultiphysics.KratosUnittest as UnitTest

class TestParallelEnvironment(UnitTest.TestCase):

    def setUp(self):
        self.world = Kratos.DataCommunicator.GetDefault()
        self.rank = self.world.Rank()
        self.size = self.world.Size()

    def tearDown(self):
        pass

    @UnitTest.skipIf(not Kratos.IsDistributedRun(), "This test is designed for distributed runs only.")
    def testDataCommunicatorRetrievalFromParallelEnvironment(self):
        default_comm = Kratos.ParallelEnvironment.GetDefaultDataCommunicator()

        self.assertTrue(Kratos.ParallelEnvironment.HasDataCommunicator("World"))

        # if we imported mpi, default should be "World" (wrapping MPI_COMM_WORLD)
        self.assertTrue(default_comm.IsDistributed())

    @UnitTest.skipIf(not Kratos.IsDistributedRun(), "This test is designed for distributed runs only.")
    def testCommunicatorCreationWithStringFromParallelEnvironment(self):
        model = Kratos.Model()
        model_part = model.CreateModelPart("MainModelPart")
        communicator = Kratos.ParallelEnvironment.CreateCommunicatorFromGlobalParallelism(model_part, "World")

        # if we imported mpi, default should be "MPICommunicator"
        self.assertTrue(communicator.GetDataCommunicator().IsDistributed())

    @UnitTest.skipIf(not Kratos.IsDistributedRun(), "This test is designed for distributed runs only.")
    def testCommunicatorCreationWithPointerFromParallelEnvironment(self):
        model = Kratos.Model()
        model_part = model.CreateModelPart("MainModelPart")
        default_comm = Kratos.ParallelEnvironment.GetDefaultDataCommunicator()
        communicator = Kratos.ParallelEnvironment.CreateCommunicatorFromGlobalParallelism(model_part, default_comm)

        # if we imported mpi, default should be "MPICommunicator"
        self.assertTrue(communicator.GetDataCommunicator().IsDistributed())

    def testDefaultRankAndSize(self):
        default_comm = Kratos.ParallelEnvironment.GetDefaultDataCommunicator()
        original_default_name = Kratos.ParallelEnvironment.GetDefaultDataCommunicatorName()
        self.assertEqual(default_comm.Rank(), Kratos.ParallelEnvironment.GetDefaultRank())
        self.assertEqual(default_comm.Size(), Kratos.ParallelEnvironment.GetDefaultSize())

        # Change the default, see if it still works
        Kratos.ParallelEnvironment.SetDefaultDataCommunicator("Serial")
        new_default = Kratos.ParallelEnvironment.GetDefaultDataCommunicator()
        self.assertEqual(new_default.Rank(), 0)
        self.assertEqual(new_default.Size(), 1)

        # Restore the original default
        Kratos.ParallelEnvironment.SetDefaultDataCommunicator(original_default_name)
        self.assertEqual(default_comm.Rank(), Kratos.ParallelEnvironment.GetDefaultRank())
        self.assertEqual(default_comm.Size(), Kratos.ParallelEnvironment.GetDefaultSize())

if __name__ == "__main__":
    UnitTest.main()