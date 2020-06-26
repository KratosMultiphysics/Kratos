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

    def markForCleanUp(self,comm_name):
        self.registered_comms.append(comm_name)

    @UnitTest.skipUnless(IsDistributedRun(), "Test is distributed.")
    def testDataCommunicatorDuplication(self):
        duplicate_comm = DataCommunicatorFactory.DuplicateAndRegister(self.default_data_communicator, "Duplicate")
        self.markForCleanUp("Duplicate") # to clean up during tearDown

        self.assertEqual(duplicate_comm.Rank(), self.default_data_communicator.Rank())
        self.assertEqual(duplicate_comm.Size(), self.default_data_communicator.Size())

    @UnitTest.skipUnless(IsDistributedRun(), "Test is distributed.")
    def testDataCommunicatorSplit(self):
        rank = self.default_data_communicator.Rank()
        size = self.default_data_communicator.Size()
        split_comm = DataCommunicatorFactory.SplitAndRegister(self.default_data_communicator, rank % 2, 0, "EvenOdd")
        self.markForCleanUp("EvenOdd") # to clean up during tearDown

        expected_rank = rank // 2
        if rank % 2 == 0:
            expected_size = math.ceil(size/2)
        else:
            expected_size = math.floor(size/2)

        self.assertEqual(split_comm.Rank(), expected_rank)
        self.assertEqual(split_comm.Size(), expected_size)

    @UnitTest.skipUnless(IsDistributedRun() and ParallelEnvironment.GetDefaultSize() > 1, "Test requires at least two ranks.")
    def testDataCommunicatorCreateFromRange(self):
        rank = self.default_data_communicator.Rank()
        size = self.default_data_communicator.Size()

        # Create a communicator using all ranks except the first
        ranks = [i for i in range(1,size)]
        range_comm = DataCommunicatorFactory.CreateFromRanksAndRegister(self.default_data_communicator, ranks, "AllExceptFirst")
        self.markForCleanUp("AllExceptFirst") # to clean up during tearDown

        if rank == 0:
            self.assertTrue(range_comm.IsNullOnThisRank())
            self.assertFalse(range_comm.IsDefinedOnThisRank())
        else:
            self.assertEqual(range_comm.Rank(), rank-1)
            self.assertEqual(range_comm.Size(), size-1)

    @UnitTest.skipUnless(IsDistributedRun() and ParallelEnvironment.GetDefaultSize() > 2, "Test requires at least three ranks.")
    def testDataCommunicatorCreateUnion(self):
        rank = self.default_data_communicator.Rank()
        size = self.default_data_communicator.Size()

        # Create a communicator using all ranks except the first
        all_except_first = DataCommunicatorFactory.CreateFromRanksAndRegister(self.default_data_communicator, [i for i in range(1,size)], "AllExceptFirst")
        self.markForCleanUp("AllExceptFirst") # to clean up during tearDown
        all_except_last = DataCommunicatorFactory.CreateFromRanksAndRegister(self.default_data_communicator, [i for i in range(0,size-1)], "AllExceptLast")
        self.markForCleanUp("AllExceptLast") # to clean up during tearDown

        # Create union communicator (should contain all ranks)
        union_comm = DataCommunicatorFactory.CreateUnionAndRegister(all_except_first, all_except_last, self.default_data_communicator, "Union")
        self.markForCleanUp("Union") # to clean up during tearDown

        self.assertFalse(union_comm.IsNullOnThisRank())
        self.assertEqual(union_comm.Rank(), rank)
        self.assertEqual(union_comm.Size(), size)

    @UnitTest.skipUnless(IsDistributedRun() and ParallelEnvironment.GetDefaultSize() > 2, "Test requires at least three ranks.")
    def testDataCommunicatorCreateIntersection(self):
        rank = self.default_data_communicator.Rank()
        size = self.default_data_communicator.Size()

        # Create a communicator using all ranks except the first
        all_except_first = DataCommunicatorFactory.CreateFromRanksAndRegister(self.default_data_communicator, [i for i in range(1,size)], "AllExceptFirst")
        self.markForCleanUp("AllExceptFirst") # to clean up during tearDown
        all_except_last = DataCommunicatorFactory.CreateFromRanksAndRegister(self.default_data_communicator, [i for i in range(0,size-1)], "AllExceptLast")
        self.markForCleanUp("AllExceptLast") # to clean up during tearDown

        intersection_comm = DataCommunicatorFactory.CreateIntersectionAndRegister(
            all_except_first, all_except_last, self.default_data_communicator, "Intersection")
        self.markForCleanUp("Intersection") # to clean up during tearDown

        if rank == 0 or rank == size - 1:
            # The first and last ranks do not participate in the intersection communicator
            self.assertTrue(intersection_comm.IsNullOnThisRank())
        else:
            self.assertEqual(intersection_comm.Rank(), rank - 1 )
            self.assertEqual(intersection_comm.Size(), size - 2 )

if __name__ == "__main__":
    UnitTest.main()
