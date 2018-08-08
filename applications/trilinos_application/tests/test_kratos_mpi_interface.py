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


if __name__ == '__main__':
    KratosUnittest.main()
