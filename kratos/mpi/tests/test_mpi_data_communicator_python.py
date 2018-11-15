import KratosMultiphysics as Kratos
import KratosMultiphysics.mpi as mpi #TODO remove once the test script accepts a command-line argument.

import KratosMultiphysics.KratosUnittest as UnitTest

class TestMPIDataCommunicatorPython(UnitTest.TestCase):

    def setUp(self):
        self.world = Kratos.ParallelEnvironment.GetInstance().GetDefaultDataCommunicator()
        self.rank = self.world.Rank()
        self.size = self.world.Size()

    def tearDown(self):
        pass

    def testDataCommunicatorRetrieval(self):
        parallel_environment = Kratos.ParallelEnvironment.GetInstance()
        default_comm = parallel_environment.GetDefaultDataCommunicator()

        self.assertTrue(parallel_environment.HasDataCommunicator("World"))

        # if we imported mpi, default should be "World" (wrapping MPI_COMM_WORLD)
        self.assertTrue(default_comm.IsDistributed())

    def testReduceListOperations(self):
        root = self.size-1
        local_int = [1, self.rank, -self.rank]
        local_double = [2.0 * i for i in local_int]

        reduced_sum_int = self.world.SumInts(local_int, root)
        reduced_min_int = self.world.MinInts(local_int, root)
        reduced_max_int = self.world.MaxInts(local_int, root)

        reduced_sum_double = self.world.SumDoubles(local_double, root)
        reduced_min_double = self.world.MinDoubles(local_double, root)
        reduced_max_double = self.world.MaxDoubles(local_double, root)

        if (self.rank == root):
            self.assertEqual(reduced_sum_int, [self.size, (self.size*(self.size-1))/2, (self.size*(self.size-1))/-2])
            self.assertEqual(reduced_min_int, [1, 0, 1-self.size])
            self.assertEqual(reduced_max_int, [1, self.size-1, 0])

            self.assertEqual(reduced_sum_double, [2.0*self.size, self.size*(self.size-1), -1.0*self.size*(self.size-1)])
            self.assertEqual(reduced_min_double, [2.0, 0.0, 2.0*(1-self.size)])
            self.assertEqual(reduced_max_double, [2.0, 2.0*(self.size-1), 0.0])

    def testReduceValueOperations(self):
        root = self.size-1

        reduced_sum_int = self.world.Sum(1, root)
        reduced_min_int = self.world.Min(self.rank, root)
        reduced_max_int = self.world.Max(self.rank, root)

        reduced_sum_double = self.world.Sum(2.0, root)
        reduced_min_double = self.world.Min(-1.0*self.rank, root)
        reduced_max_double = self.world.Max(-1.0*self.rank, root)

        if (self.rank == root):
            self.assertEqual(reduced_sum_int, self.size)
            self.assertEqual(reduced_min_int, 0)
            self.assertEqual(reduced_max_int, self.size-1)

            self.assertEqual(reduced_sum_double, 2.0*self.size)
            self.assertEqual(reduced_min_double, -1.0*(self.size-1))
            self.assertEqual(reduced_max_double, 0.0)

    def testAllReduceListOperations(self):
        local_int = [1, self.rank, -self.rank]
        local_double = [2.0 * i for i in local_int]

        reduced_sum_int = self.world.SumAllInts(local_int)
        reduced_min_int = self.world.MinAllInts(local_int)
        reduced_max_int = self.world.MaxAllInts(local_int)

        reduced_sum_double = self.world.SumAllDoubles(local_double)
        reduced_min_double = self.world.MinAllDoubles(local_double)
        reduced_max_double = self.world.MaxAllDoubles(local_double)

        self.assertEqual(reduced_sum_int, [self.size, (self.size*(self.size-1))/2, (self.size*(self.size-1))/-2])
        self.assertEqual(reduced_min_int, [1, 0, 1-self.size])
        self.assertEqual(reduced_max_int, [1, self.size-1, 0])

        self.assertEqual(reduced_sum_double, [2.0*self.size, self.size*(self.size-1), -1.0*self.size*(self.size-1)])
        self.assertEqual(reduced_min_double, [2.0, 0.0, 2.0*(1-self.size)])
        self.assertEqual(reduced_max_double, [2.0, 2.0*(self.size-1), 0.0])

    def testAllReduceValueOperations(self):
        reduced_sum_int = self.world.SumAll(1)
        reduced_min_int = self.world.MinAll(self.rank)
        reduced_max_int = self.world.MaxAll(self.rank)

        reduced_sum_double = self.world.SumAll(2.0)
        reduced_min_double = self.world.MinAll(-1.0*self.rank)
        reduced_max_double = self.world.MaxAll(-1.0*self.rank)

        self.assertEqual(reduced_sum_int, self.size)
        self.assertEqual(reduced_min_int, 0)
        self.assertEqual(reduced_max_int, self.size-1)

        self.assertEqual(reduced_sum_double, 2.0*self.size)
        self.assertEqual(reduced_min_double, -1.0*(self.size-1))
        self.assertEqual(reduced_max_double, 0.0)

    def testScanSumOperations(self):
        partial_sum_int = self.world.ScanSum(-1)
        partial_sum_double = self.world.ScanSum(2.0)
        partial_sum_int_list = self.world.ScanSumInts([1,-1])
        partial_sum_double_list = self.world.ScanSumDoubles([2.0,-2.0])

        self.assertEqual(partial_sum_int, -self.rank-1)
        self.assertEqual(partial_sum_double, 2.0*(self.rank+1))
        self.assertEqual(partial_sum_int_list, [self.rank+1, -self.rank-1])
        self.assertEqual(partial_sum_double_list, [2.0*(self.rank+1), -2.0*(self.rank+1)])

    @UnitTest.skipIf(Kratos.ParallelEnvironment.GetInstance().GetDefaultDataCommunicator().Size() < 2, "This test does not work for a single process.")
    def testSendRecvOperations(self):

        if self.rank + 1 != self.size:
            send_rank = self.rank + 1
        else:
            send_rank = 0

        if self.rank != 0:
            recv_rank = self.rank - 1
        else:
            recv_rank = self.size - 1

        send_ints = [self.rank, self.rank]
        send_doubles = [2.0*self.rank, 2.0*self.rank]

        recv_ints    = self.world.SendRecvInts(   send_ints,    send_rank, recv_rank)
        recv_doubles = self.world.SendRecvDoubles(send_doubles, send_rank, recv_rank)

        if self.rank == 0:
            self.assertEqual(recv_ints, [self.size-1, self.size-1])
            self.assertEqual(recv_doubles, [2.0*self.size-2, 2.0*self.size-2])
        else:
            self.assertEqual(recv_ints, [self.rank-1, self.rank-1])
            self.assertEqual(recv_doubles, [2.0*self.rank-2, 2.0*self.rank-2])

    def testBroadcastOperations(self):
        source_rank = 0
        broadcast_int = self.world.Broadcast(self.rank, source_rank)
        broadcast_double = self.world.Broadcast(2.0*self.rank, source_rank)
        broadcast_int_list = self.world.BroadcastInts([self.rank,-self.rank], source_rank)
        broadcast_double_list = self.world.BroadcastDoubles([2.0*self.rank,-2.0*self.rank], source_rank)

        self.assertEqual(broadcast_int, source_rank)
        self.assertEqual(broadcast_double, 2.0*source_rank)
        self.assertEqual(broadcast_int_list, [source_rank, -source_rank])
        self.assertEqual(broadcast_double_list, [2.0*source_rank, -2.0*source_rank])

    def testScatterOperations(self):
        source_rank = self.size-1

        if self.rank == source_rank:
            scatter_source_ints = [i for i in range(self.size)]
            scatter_source_double = [i for i in range(2*self.size)]
        else:
            scatter_source_ints = list()
            scatter_source_double = list()

        scatter_int_list = self.world.ScatterInts(scatter_source_ints, source_rank)
        scatter_double_list = self.world.ScatterDoubles(scatter_source_double, source_rank)

        self.assertEqual(scatter_int_list, [self.rank,])
        self.assertEqual(scatter_double_list, [2.0*self.rank,2.0*self.rank+1])

    def testScattervOperations(self):
        source_rank = self.size-1

        if self.rank == source_rank:
            scatter_source_ints = [ [i, i+1, i+2] for i in range(self.size)]
            scatter_source_double = [ [1.0 + i, 2.0 + i , 3.0 + i] for i in range(self.size)]
        else:
            scatter_source_ints = list()
            scatter_source_double = list()

        scatter_int_list = self.world.ScattervInts(scatter_source_ints, source_rank)
        scatter_double_list = self.world.ScattervDoubles(scatter_source_double, source_rank)

        self.assertEqual(scatter_int_list, [self.rank,self.rank+1, self.rank+2])
        self.assertEqual(scatter_double_list, [1.0+self.rank, 2.0+self.rank, 3.0+self.rank])

    def testGatherOperations(self):
        destination_rank = self.size-1

        message_int = [self.rank,]
        message_double = [2.0*self.rank,]

        gathered_ints = self.world.GatherInts(message_int, destination_rank)
        gathered_doubles = self.world.GatherDoubles(message_double, destination_rank)

        if self.rank == destination_rank:
            self.assertEqual(gathered_ints, [i for i in range(self.size)])
            self.assertEqual(gathered_doubles, [2.0*i for i in range(self.size)])
        else:
            self.assertEqual(gathered_ints, [])
            self.assertEqual(gathered_doubles, [])

    def testGathervOperations(self):
        destination_rank = 0

        message_int = [self.rank, self.rank+1, self.rank+2]
        message_double = [1.0 + self.rank, 2.0 + self.rank, 3.0 + self.rank]

        gathered_ints = self.world.GathervInts(message_int, destination_rank)
        gathered_doubles = self.world.GathervDoubles(message_double, destination_rank)

        if self.rank == destination_rank:
            self.assertEqual(gathered_ints, [ [i, i+1, i+2] for i in range(self.size)])
            self.assertEqual(gathered_doubles, [ [1.0 + i, 2.0 + i , 3.0 + i] for i in range(self.size)])
        else:
            self.assertEqual(gathered_ints, [[] for i in range(self.size)])
            self.assertEqual(gathered_doubles, [[] for i in range(self.size)])

    def testAllGatherOperations(self):
        message_int = [self.rank,]
        message_double = [2.0*self.rank,]

        gathered_ints = self.world.AllGatherInts(message_int)
        gathered_doubles = self.world.AllGatherDoubles(message_double)

        self.assertEqual(gathered_ints, [i for i in range(self.size)])
        self.assertEqual(gathered_doubles, [2.0*i for i in range(self.size)])


if __name__ == "__main__":
    UnitTest.main()