import KratosMultiphysics as Kratos

import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestMPIDataCommunicatorPython(KratosUnittest.TestCase):

    def setUp(self):
        self.world: Kratos.DataCommunicator = Kratos.Testing.GetDefaultDataCommunicator()
        self.rank = self.world.Rank()
        self.size = self.world.Size()

    def tearDown(self):
        pass

    @KratosUnittest.skipIf(not Kratos.IsDistributedRun(), "This test is designed for distributed runs only.")
    def testDataCommunicatorRetrievalFromDataCommunicator(self):
        default_comm = Kratos.Testing.GetDefaultDataCommunicator()

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

        if self.rank == root:
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

        if self.rank == root:
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

    def testAllReduceLocValueOperations(self):
        reduced_min_int = self.world.MinLocAll(self.rank)
        reduced_max_int = self.world.MaxLocAll(self.rank)

        reduced_min_double = self.world.MinLocAll(-1.0*self.rank)
        reduced_max_double = self.world.MaxLocAll(-1.0*self.rank)

        self.assertEqual(reduced_min_int[0], 0)
        self.assertEqual(reduced_min_int[1], 0)
        self.assertEqual(reduced_max_int[0], self.size-1)
        self.assertEqual(reduced_max_int[1], self.size-1)

        self.assertEqual(reduced_min_double[0], -1.0*(self.size-1))
        self.assertEqual(reduced_min_double[1], self.size-1)
        self.assertEqual(reduced_max_double[0], 0.0)
        self.assertEqual(reduced_max_double[1], 0)

    def testScanSumOperations(self):
        partial_sum_int = self.world.ScanSum(-1)
        partial_sum_double = self.world.ScanSum(2.0)
        partial_sum_int_list = self.world.ScanSumInts([1,-1])
        partial_sum_double_list = self.world.ScanSumDoubles([2.0,-2.0])

        self.assertEqual(partial_sum_int, -self.rank-1)
        self.assertEqual(partial_sum_double, 2.0*(self.rank+1))
        self.assertEqual(partial_sum_int_list, [self.rank+1, -self.rank-1])
        self.assertEqual(partial_sum_double_list, [2.0*(self.rank+1), -2.0*(self.rank+1)])

    @KratosUnittest.skipIf(Kratos.ParallelEnvironment.GetDefaultDataCommunicator().Size() < 2, "This test does not work for a single process.")
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

    @KratosUnittest.skipIf(Kratos.ParallelEnvironment.GetDefaultDataCommunicator().Size() < 2, "This test does not work for a single process.")
    def testSendRecvString(self):

        if self.rank + 1 != self.size:
            send_rank = self.rank + 1
        else:
            send_rank = 0

        if self.rank != 0:
            recv_rank = self.rank - 1
        else:
            recv_rank = self.size - 1

        send_string = "Hello from rank {0}".format(self.rank)
        send_lengths = [len(send_string),]

        recv_lengths = self.world.SendRecvInts(   send_lengths, send_rank, recv_rank)
        recv_string  = self.world.SendRecvString( send_string,  send_rank, recv_rank)

        self.assertEqual(len(recv_string), recv_lengths[0])
        self.assertEqual(recv_string, "Hello from rank {0}".format(recv_rank))

    def testBroadcastOperations(self):
        source_rank = 0
        broadcast_int = self.world.Broadcast(self.rank, source_rank)
        broadcast_double = self.world.Broadcast(2.0*self.rank, source_rank)
        broadcast_string = self.world.Broadcast(str(self.rank), source_rank)
        broadcast_vector = self.world.Broadcast(Kratos.Vector(3, self.rank + 1), source_rank)
        broadcast_matrix = self.world.Broadcast(Kratos.Matrix(3, 3, self.rank + 1), source_rank)
        broadcast_int_list = self.world.BroadcastInts([self.rank,-self.rank], source_rank)
        broadcast_double_list = self.world.BroadcastDoubles([2.0*self.rank,-2.0*self.rank], source_rank)
        broadcast_string_list = self.world.BroadcastStrings([f"{self.rank}_{i}" for i in range(3)], source_rank)

        self.assertEqual(broadcast_int, source_rank)
        self.assertEqual(broadcast_double, 2.0*source_rank)
        self.assertEqual(broadcast_string, str(source_rank))
        self.assertEqual(broadcast_int_list, [source_rank, -source_rank])
        self.assertEqual(broadcast_double_list, [2.0*source_rank, -2.0*source_rank])
        self.assertEqual(broadcast_string_list, [f"{source_rank}_{i}" for i in range(3)])
        self.assertVectorAlmostEqual(broadcast_vector, Kratos.Vector(3, source_rank + 1), 12)
        self.assertMatrixAlmostEqual(broadcast_matrix, Kratos.Matrix(3, 3, source_rank + 1), 12)

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

    def testAllGathervOperations(self):
        message_int = [self.rank, self.rank+1, self.rank+2]
        message_double = [1.0 + self.rank, 2.0 + self.rank, 3.0 + self.rank]

        gathered_ints = self.world.AllGathervInts(message_int)
        gathered_doubles = self.world.AllGathervDoubles(message_double)

        self.assertEqual(gathered_ints, [ [i, i+1, i+2] for i in range(self.size)])
        self.assertEqual(gathered_doubles, [ [1.0 + i, 2.0 + i , 3.0 + i] for i in range(self.size)])

    def test_CastingTypes(self):
        n = self.world.Size()

        def check_value_in_rank(ref_value, value, ranks_to_check: 'list[int]'):
            self.assertTrue(isinstance(value, type(ref_value)))
            if self.rank in ranks_to_check:
                if isinstance(ref_value, int) or isinstance(ref_value, float):
                    self.assertEqual(ref_value, value)
                elif isinstance(ref_value, Kratos.Array3) or isinstance(ref_value, Kratos.Array4) or isinstance(ref_value, Kratos.Array6) or isinstance(ref_value, Kratos.Array9) or isinstance(ref_value, Kratos.Vector):
                    self.assertVectorAlmostEqual(ref_value, value)
                elif isinstance(ref_value, Kratos.Matrix):
                    self.assertMatrixAlmostEqual(ref_value, value)
                elif isinstance(ref_value, list):
                    check_value_in_rank(ref_value[0], value[0], ranks_to_check)

        def simple_reduce_check(method, rank: int, ref_value: float, ranks_to_check: 'list[int]'):
            check_value_in_rank(method(self.rank+1, rank), int(ref_value), ranks_to_check)
            check_value_in_rank(method(float(self.rank+1), rank), float(ref_value), ranks_to_check)
            check_value_in_rank(method(Kratos.Array3(self.rank+1), rank), Kratos.Array3(ref_value), ranks_to_check)
            check_value_in_rank(method(Kratos.Array4(self.rank+1), rank), Kratos.Array4(ref_value), ranks_to_check)
            check_value_in_rank(method(Kratos.Array6(self.rank+1), rank), Kratos.Array6(ref_value), ranks_to_check)
            check_value_in_rank(method(Kratos.Array9(self.rank+1), rank), Kratos.Array9(ref_value), ranks_to_check)
            check_value_in_rank(method(Kratos.Vector(2, self.rank+1), rank), Kratos.Vector(2, ref_value), ranks_to_check)
            check_value_in_rank(method(Kratos.Matrix(2, 2, self.rank+1), rank), Kratos.Matrix(2, 2, ref_value), ranks_to_check)

            check_value_in_rank(getattr(self.world, method.__name__ + "Ints")([self.rank+1], rank), [int(ref_value)], ranks_to_check)
            check_value_in_rank(getattr(self.world, method.__name__ + "Doubles")([float(self.rank+1)], rank), [float(ref_value)], ranks_to_check)
            check_value_in_rank(getattr(self.world, method.__name__ + "Array3s")([Kratos.Array3(self.rank+1)], rank), [Kratos.Array3(ref_value)], ranks_to_check)
            check_value_in_rank(getattr(self.world, method.__name__ + "Array4s")([Kratos.Array4(self.rank+1)], rank), [Kratos.Array4(ref_value)], ranks_to_check)
            check_value_in_rank(getattr(self.world, method.__name__ + "Array6s")([Kratos.Array6(self.rank+1)], rank), [Kratos.Array6(ref_value)], ranks_to_check)
            check_value_in_rank(getattr(self.world, method.__name__ + "Array9s")([Kratos.Array9(self.rank+1)], rank), [Kratos.Array9(ref_value)], ranks_to_check)
            check_value_in_rank(getattr(self.world, method.__name__ + "Vectors")([Kratos.Vector(2, self.rank+1)], rank), [Kratos.Vector(2, ref_value)], ranks_to_check)
            check_value_in_rank(getattr(self.world, method.__name__ + "Matrices")([Kratos.Matrix(2, 2, self.rank+1)], rank), [Kratos.Matrix(2, 2, ref_value)], ranks_to_check)

        def simple_all_reduce_check(method, ref_value: float):
            ranks_to_check = [i for i in range(n)]
            check_value_in_rank(method(self.rank+1), int(ref_value), ranks_to_check)
            check_value_in_rank(method(float(self.rank+1)), float(ref_value), ranks_to_check)
            check_value_in_rank(method(Kratos.Array3(self.rank+1)), Kratos.Array3(ref_value), ranks_to_check)
            check_value_in_rank(method(Kratos.Array3(self.rank+1)), Kratos.Array3(ref_value), ranks_to_check)
            check_value_in_rank(method(Kratos.Array3(self.rank+1)), Kratos.Array3(ref_value), ranks_to_check)
            check_value_in_rank(method(Kratos.Array3(self.rank+1)), Kratos.Array3(ref_value), ranks_to_check)
            check_value_in_rank(method(Kratos.Vector(2, self.rank+1)), Kratos.Vector(2, ref_value), ranks_to_check)
            check_value_in_rank(method(Kratos.Matrix(2, 2, self.rank+1)), Kratos.Matrix(2, 2, ref_value), ranks_to_check)

            check_value_in_rank(getattr(self.world, method.__name__ + "Ints")([self.rank+1]), [int(ref_value)], ranks_to_check)
            check_value_in_rank(getattr(self.world, method.__name__ + "Doubles")([float(self.rank+1)]), [float(ref_value)], ranks_to_check)
            check_value_in_rank(getattr(self.world, method.__name__ + "Array3s")([Kratos.Array3(self.rank+1)]), [Kratos.Array3(ref_value)], ranks_to_check)
            check_value_in_rank(getattr(self.world, method.__name__ + "Array4s")([Kratos.Array4(self.rank+1)]), [Kratos.Array4(ref_value)], ranks_to_check)
            check_value_in_rank(getattr(self.world, method.__name__ + "Array6s")([Kratos.Array6(self.rank+1)]), [Kratos.Array6(ref_value)], ranks_to_check)
            check_value_in_rank(getattr(self.world, method.__name__ + "Array9s")([Kratos.Array9(self.rank+1)]), [Kratos.Array9(ref_value)], ranks_to_check)
            check_value_in_rank(getattr(self.world, method.__name__ + "Vectors")([Kratos.Vector(2, self.rank+1)]), [Kratos.Vector(2, ref_value)], ranks_to_check)
            check_value_in_rank(getattr(self.world, method.__name__ + "Matrices")([Kratos.Matrix(2, 2, self.rank+1)]), [Kratos.Matrix(2, 2, ref_value)], ranks_to_check)

        simple_reduce_check(self.world.Sum, 0, n*(n+1)/2, [0])
        simple_reduce_check(self.world.Min, 0, 1, [0])
        simple_reduce_check(self.world.Max, 0, n, [0])
        simple_reduce_check(self.world.Broadcast, 0, 1, [i for i in range(n)])
        simple_all_reduce_check(self.world.SumAll, n*(n+1)/2)
        simple_all_reduce_check(self.world.MinAll, 1)
        simple_all_reduce_check(self.world.MaxAll, n)

if __name__ == "__main__":
    Kratos.Logger.GetDefaultOutput().SetSeverity(Kratos.Logger.Severity.WARNING)
    KratosUnittest.main()
