from __future__ import print_function, absolute_import, division

import KratosMultiphysics
import KratosMultiphysics.mpi as mpi

class LegacyMPICommInterface(object):

    def __init__(self, comm):
        self.comm = comm

    def barrier(self):
        self.comm.Barrier()

class LegacyMPIPythonInterface(object):

    def __init__(self):
        self.__world = KratosMultiphysics.ParallelEnvironment.GetDataCommunicator("World")
        self.world = LegacyMPICommInterface(self.__world)
        self.rank = self.__world.Rank()
        self.size = self.__world.Size()

    def broadcast_double(self, comm_wrapper, value, source_rank):
        return comm_wrapper.comm.Broadcast(value, source_rank)

    def broadcast_int(self, comm_wrapper, value, source_rank):
        return comm_wrapper.comm.Broadcast(value, source_rank)

    def max_double(self, comm_wrapper, value, reduce_rank):
        return comm_wrapper.comm.Max(value, reduce_rank)

    def max_int(self, comm_wrapper, value, reduce_rank):
        return comm_wrapper.comm.Max(value, reduce_rank)

    def min_double(self, comm_wrapper, value, reduce_rank):
        return comm_wrapper.comm.Min(value, reduce_rank)

    def min_int(self, comm_wrapper, value, reduce_rank):
        return comm_wrapper.comm.Min(value, reduce_rank)

    def sum_double(self, comm_wrapper, value, reduce_rank):
        return comm_wrapper.comm.Sum(value, reduce_rank)

    def sum_int(self, comm_wrapper, value, reduce_rank):
        return comm_wrapper.comm.Sum(value, reduce_rank)

    def max_all_double(self, comm_wrapper, value):
        return comm_wrapper.comm.MaxAll(value)

    def max_all_int(self, comm_wrapper, value):
        return comm_wrapper.comm.MaxAll(value)

    def min_all_double(self, comm_wrapper, value):
        return comm_wrapper.comm.MinAll(value)

    def min_all_int(self, comm_wrapper, value):
        return comm_wrapper.comm.MinAll(value)

    def sum_all_double(self, comm_wrapper, value):
        return comm_wrapper.comm.SumAll(value)

    def sum_all_int(self, comm_wrapper, value):
        return comm_wrapper.comm.SumAll(value)

    def scatter_double(self, comm_wrapper, local_values, source_rank):
        output_list = comm_wrapper.comm.ScatterDoubles(local_values, source_rank)
        return output_list[0]

    def scatter_int(self, comm_wrapper, local_values, source_rank):
        output_list = comm_wrapper.comm.ScatterInts(local_values, source_rank)
        return output_list[0]

    def scatterv_double(self, comm_wrapper, local_values, source_rank):
        return comm_wrapper.comm.ScattervDoubles(local_values, source_rank)

    def scatterv_int(self, comm_wrapper, local_values, source_rank):
        return comm_wrapper.comm.ScattervInts(local_values, source_rank)

    def gather_double(self, comm_wrapper, local_value, gather_rank):
        wraper_list = [local_value,]
        return comm_wrapper.comm.GatherDoubles(wraper_list, gather_rank)

    def gather_int(self, comm_wrapper, local_value, gather_rank):
        wraper_list = [local_value,]
        return comm_wrapper.comm.GatherInts(wraper_list, gather_rank)

    def gatherv_double(self, comm_wrapper, local_values, gather_rank):
        output = comm_wrapper.comm.GathervDoubles(local_values, gather_rank)
        if self.__world.Rank() == gather_rank:
            return output
        else:
            return []

    def gatherv_int(self, comm_wrapper, local_values, gather_rank):
        output = comm_wrapper.comm.GathervInts(local_values, gather_rank)
        if self.__world.Rank() == gather_rank:
            return output
        else:
            return []

    def allgather_double(self, comm_wrapper, local_value):
        wraper_list = [local_value,]
        return comm_wrapper.comm.AllGatherDoubles(wraper_list)

    def allgather_int(self, comm_wrapper, local_value):
        wraper_list = [local_value,]
        return comm_wrapper.comm.AllGatherInts(wraper_list)


