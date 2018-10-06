from __future__ import print_function, absolute_import, division
import KratosMultiphysics
from KratosMultiphysics.mpi import mpi
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestKratosMPIInterface(KratosUnittest.TestCase):
    def test_gather_int(self):
        my_rank = mpi.rank
        rank_to_gather_on = 0

        my_value = my_rank # this is an int

        gathered_values = mpi.gather(mpi.world, my_value, rank_to_gather_on)

        self.assertEqual(type(gathered_values), list)

        if my_rank == rank_to_gather_on:
            comm_size = mpi.size
            self.assertEqual(len(gathered_values), comm_size)
            exp_value = int((comm_size-1)*(comm_size)/2) # n(n+1)/2
            self.assertEqual(sum(gathered_values), exp_value)
        else:
            self.assertEqual(gathered_values, [])

    def test_gather_double(self):
        my_rank = mpi.rank
        rank_to_gather_on = 0

        my_value = my_rank+0.5 # this is a double

        gathered_values = mpi.gather(mpi.world, my_value, rank_to_gather_on)

        self.assertEqual(type(gathered_values), list)

        if my_rank == rank_to_gather_on:
            comm_size = mpi.size
            self.assertEqual(len(gathered_values), comm_size)
            exp_value = (comm_size-1)*(comm_size)/2 + comm_size*0.5 # n(n+1)/2 + comm_size*0.5
            self.assertEqual(sum(gathered_values), exp_value)
        else:
            self.assertEqual(gathered_values, [])

    def test_gatherv(self):
        my_rank = mpi.rank
        rank_to_gather_on = 0

        local_list = []

        for i in range(my_rank+1):
            local_list.append(i+1.3)

        gathered_lists = mpi.gatherv(mpi.world, local_list, rank_to_gather_on) # this is a list of lists, one from each rank

        self.assertEqual(type(gathered_lists), list)
        for sublist in gathered_lists:
            self.assertEqual(type(sublist), list)

        if my_rank == rank_to_gather_on:
            comm_size = mpi.size
            self.assertEqual(len(gathered_lists), comm_size) # checking the number of lists
            for i in range(comm_size):
                self.assertEqual(len(gathered_lists[i]), i+1) # checking the length of the lists returned from the ranks
                n = i+2
                exp_value = (n-1)*(n)/2 +(i+1)*0.3 # n(n+1)/2
                self.assertAlmostEqual(sum(gathered_lists[i]), exp_value)
        else:
            self.assertEqual(gathered_lists, [])

    def test_allgather_int(self):
        my_value = mpi.rank # this is an int

        gathered_values = mpi.allgather(mpi.world, my_value)

        self.assertEqual(type(gathered_values), list)

        comm_size = mpi.size
        self.assertEqual(len(gathered_values), comm_size)
        exp_value = (comm_size-1)*(comm_size)/2 # n(n+1)/2
        self.assertEqual(sum(gathered_values), exp_value)

    def test_allgather_double(self):
        my_value = mpi.rank+0.5 # this is a double

        gathered_values = mpi.allgather(mpi.world, my_value)

        self.assertEqual(type(gathered_values), list)

        comm_size = mpi.size
        self.assertEqual(len(gathered_values), comm_size)
        exp_value = (comm_size-1)*(comm_size)/2 + comm_size*0.5 # n(n+1)/2 + comm_size*0.5
        self.assertEqual(sum(gathered_values), exp_value)

    def test_scatter_int(self):
        rank_to_scatter_from = 0
        my_rank = mpi.rank
        offset = 21
        if my_rank == rank_to_scatter_from:
            vals_to_scatter = list(range(offset, mpi.size+offset))
        else:
            vals_to_scatter = list()
        scattered_val = mpi.scatter_int(mpi.world, vals_to_scatter, rank_to_scatter_from)
        self.assertEqual(scattered_val, my_rank+offset)

    def test_scatter_double(self):
        rank_to_scatter_from = 0
        my_rank = mpi.rank
        size = mpi.size
        offset_1 = 21.4
        offset_2 = 0.22
        vals_to_scatter = list()
        if my_rank == rank_to_scatter_from:
            for i in range(size):
                vals_to_scatter.append(offset_1 + offset_2 + i)
        scattered_val = mpi.scatter_double(mpi.world, vals_to_scatter, rank_to_scatter_from)
        self.assertAlmostEqual(scattered_val, mpi.rank+offset_1+offset_2)

    def test_scatterv_int(self):
        rank_to_scatter_from = 0
        my_rank = mpi.rank
        min_list_size = 3
        if my_rank == rank_to_scatter_from:
            comm_size = mpi.size
            vals_to_scatter = list()
            for i in range(comm_size):
                rank_list = list(range(i+min_list_size))
                vals_to_scatter.append(rank_list)
        else:
            vals_to_scatter = list(list())
        scattered_list = mpi.scatterv_int(mpi.world, vals_to_scatter, rank_to_scatter_from)

        self.assertEqual(type(scattered_list), list)
        self.assertEqual(len(scattered_list), my_rank+min_list_size)

        my_rank_list = list(range(my_rank+min_list_size))

        for recv_val, exp_val in zip(scattered_list, my_rank_list):
            self.assertEqual(recv_val, exp_val)

    def test_scatterv_double(self):
        rank_to_scatter_from = 0
        my_rank = mpi.rank
        min_list_size = 3
        if my_rank == rank_to_scatter_from:
            comm_size = mpi.size
            vals_to_scatter = list()
            for i in range(comm_size):
                range_list = list(range(i+min_list_size))
                rank_list = DoubleScatterVList(range_list)
                vals_to_scatter.append(rank_list)
        else:
            vals_to_scatter = list(list())
        scattered_list = mpi.scatterv_double(mpi.world, vals_to_scatter, rank_to_scatter_from)

        self.assertEqual(type(scattered_list), list)
        self.assertEqual(len(scattered_list), my_rank+min_list_size, msg="Rank: "+str(my_rank) + " ; " + str(scattered_list))

        my_range_list = list(range(my_rank+min_list_size))
        my_rank_list = DoubleScatterVList(my_range_list)

        for recv_val, exp_val in zip(scattered_list, my_rank_list):
            self.assertAlmostEqual(recv_val, exp_val)

    def test_broadcast_int(self):
        val_to_broadcast = (mpi.rank+1)*5
        rank_to_broadcast_from = 0

        broadcasted_value = mpi.broadcast(mpi.world, val_to_broadcast, rank_to_broadcast_from)

        exp_val = (rank_to_broadcast_from+1)*5

        self.assertEqual(exp_val, broadcasted_value)

    def test_broadcast_double(self):
        val_to_broadcast = (mpi.rank+1)*5.333
        rank_to_broadcast_from = 0

        broadcasted_value = mpi.broadcast(mpi.world, val_to_broadcast, rank_to_broadcast_from)

        exp_val = (rank_to_broadcast_from+1)*5.333

        self.assertAlmostEqual(exp_val, broadcasted_value)

    def test_reduce_int_max(self):
        self._execute_reduction_test_int_max(is_allreduce=False)

    def test_allreduce_int_max(self):
        self._execute_reduction_test_int_max(is_allreduce=True)

    def test_reduce_int_min(self):
        self._execute_reduction_test_int_min(is_allreduce=False)

    def test_allreduce_int_min(self):
        self._execute_reduction_test_int_min(is_allreduce=True)

    def test_reduce_int_sum(self):
        self._execute_reduction_test_int_sum(is_allreduce=False)

    def test_allreduce_int_sum(self):
        self._execute_reduction_test_int_sum(is_allreduce=True)

    def _execute_reduction_test_int_max(self, is_allreduce):
        my_rank = mpi.rank
        rank_to_reduce_on = 0

        my_value = my_rank + 2 # this is an int

        if is_allreduce:
            max_val = mpi.allreduce(mpi.world, my_value, mpi.MPI_op.MAX)
        else:
            max_val = mpi.reduce(mpi.world, my_value, rank_to_reduce_on, mpi.MPI_op.MAX)

        if is_allreduce or my_rank == rank_to_reduce_on:
            self.assertEqual(mpi.size+1, max_val)

    def _execute_reduction_test_int_min(self, is_allreduce):
        my_rank = mpi.rank
        rank_to_reduce_on = 0

        my_value = my_rank - 2 # this is an int

        if is_allreduce:
            min_val = mpi.allreduce(mpi.world, my_value, mpi.MPI_op.MIN)
        else:
            min_val = mpi.reduce(mpi.world, my_value, rank_to_reduce_on, mpi.MPI_op.MIN)

        if is_allreduce or my_rank == rank_to_reduce_on:
            self.assertEqual(-2, min_val)

    def _execute_reduction_test_int_sum(self, is_allreduce):
        my_rank = mpi.rank
        rank_to_reduce_on = 0

        my_value = my_rank+1 # this is an int

        if is_allreduce:
            sum_val = mpi.allreduce(mpi.world, my_value, mpi.MPI_op.SUM)
        else:
            sum_val = mpi.reduce(mpi.world, my_value, rank_to_reduce_on, mpi.MPI_op.SUM)

        if is_allreduce or my_rank == rank_to_reduce_on:
            comm_size = mpi.size
            exp_val = int((comm_size)*(comm_size+1)/2)
            self.assertEqual(exp_val, sum_val)

def DoubleScatterVList(range_list):
    return [(2.125**x)/(2**(x-6)) for x in range_list]


if __name__ == '__main__':
    KratosUnittest.main()
