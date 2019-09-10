from KratosMultiphysics import ParallelEnvironment, IsDistributedRun
if IsDistributedRun():
    from KratosMultiphysics.mpi import DataCommunicatorFactory
import KratosMultiphysics.KratosUnittest as UnitTest

import math

class TestDataCommunicatorFactory(UnitTest.TestCase):

    def setUp(self):
        self.registered_comms = []
        self.default_data_communicator = ParallelEnvironment.GetDefaultDataCommunicator()
        self.original_default = ParallelEnvironment.GetDefaultDataCommunicatorName()

    def tearDown(self):
        if len(self.registered_comms) > 0:
            ParallelEnvironment.SetDefaultDataCommunicator(self.original_default)
            for comm_name in self.registered_comms:
                ParallelEnvironment.UnregisterDataCommunicator(comm_name)

    @UnitTest.skipUnless(IsDistributedRun(), "Test is distributed.")
    def testDataCommunicatorDuplication(self):
        duplicate_comm = DataCommunicatorFactory.DuplicateAndRegister(self.default_data_communicator, "Duplicate")
        self.registered_comms.append("Duplicate") # to clean up during tearDown

        self.assertEqual(duplicate_comm.Rank(), self.default_data_communicator.Rank())
        self.assertEqual(duplicate_comm.Size(), self.default_data_communicator.Size())

    @UnitTest.skipUnless(IsDistributedRun(), "Test is distributed.")
    def testDataCommunicatorSplit(self):
        rank = self.default_data_communicator.Rank()
        size = self.default_data_communicator.Size()
        split_comm = DataCommunicatorFactory.SplitAndRegister(self.default_data_communicator, rank % 2, 0, "EvenOdd")
        self.registered_comms.append("EvenOdd") # to clean up during tearDown

        expected_rank = rank // 2
        if rank % 2 == 0:
            expected_size = math.ceil(size/2)
        else:
            expected_size = math.floor(size/2)

        self.assertEqual(split_comm.Rank(), expected_rank)
        self.assertEqual(split_comm.Size(), expected_size)

if __name__ == "__main__":
    UnitTest.main()
