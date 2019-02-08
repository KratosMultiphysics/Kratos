from __future__ import print_function, absolute_import, division

import KratosMultiphysics

class LegacyMPICommInterface(object):

    def __init__(self, comm):
        self.comm = comm

    def barrier(self):
        self.GetWithDeprecationWarning("barrier").Barrier()

    def GetWithDeprecationWarning(self, func_name):
        if self.comm.Rank() == 0:
            msg = [
                "Calling the legacy MPI Python interface function ",
                func_name,
                ". Please use a DataCommunicator instance instead.\n",
                "DataCommunicators can be obtained calling\n",
                "  KratosMultiphysics.ParallelEnvironment.GetDataCommunicator(\"World\")\n",
                "(for the DataCommunicator based on MPI_COMM_WORLD) or \n",
                "  KratosMultiphysics.DataCommunicator.GetDefault()\n",
                "See the tests in kratos/mpi/tests/test_mpi_data_communicator_python.py for examples of usage."
            ]
            KratosMultiphysics.Logger.PrintWarning("MPI Interface","".join(msg))
        return self.comm


class LegacyMPIPythonInterface(object):

    def __init__(self):
        self.__world = KratosMultiphysics.ParallelEnvironment.GetDataCommunicator("World")
        self.world = LegacyMPICommInterface(self.__world)
        self.rank = self.__world.Rank()
        self.size = self.__world.Size()

    def broadcast_double(self, comm_wrapper, value, source_rank):
        return comm_wrapper.GetWithDeprecationWarning("broadcast_double").Broadcast(value, source_rank)

    def broadcast_int(self, comm_wrapper, value, source_rank):
        return comm_wrapper.GetWithDeprecationWarning("broadcast_int").Broadcast(value, source_rank)

    def max_double(self, comm_wrapper, value, reduce_rank):
        return comm_wrapper.GetWithDeprecationWarning("max_double").Max(value, reduce_rank)

    def max_int(self, comm_wrapper, value, reduce_rank):
        return comm_wrapper.GetWithDeprecationWarning("max_int").Max(value, reduce_rank)

    def min_double(self, comm_wrapper, value, reduce_rank):
        return comm_wrapper.GetWithDeprecationWarning("min_double").Min(value, reduce_rank)

    def min_int(self, comm_wrapper, value, reduce_rank):
        return comm_wrapper.GetWithDeprecationWarning("min_int").Min(value, reduce_rank)

    def sum_double(self, comm_wrapper, value, reduce_rank):
        return comm_wrapper.GetWithDeprecationWarning("sum_double").Sum(value, reduce_rank)

    def sum_int(self, comm_wrapper, value, reduce_rank):
        return comm_wrapper.GetWithDeprecationWarning("sum_int").Sum(value, reduce_rank)

    def max_all_double(self, comm_wrapper, value):
        return comm_wrapper.GetWithDeprecationWarning("max_all_double").MaxAll(value)

    def max_all_int(self, comm_wrapper, value):
        return comm_wrapper.GetWithDeprecationWarning("max_all_int").MaxAll(value)

    def min_all_double(self, comm_wrapper, value):
        return comm_wrapper.GetWithDeprecationWarning("min_all_double").MinAll(value)

    def min_all_int(self, comm_wrapper, value):
        return comm_wrapper.GetWithDeprecationWarning("min_all_int").MinAll(value)

    def sum_all_double(self, comm_wrapper, value):
        return comm_wrapper.GetWithDeprecationWarning("sum_all_double").SumAll(value)

    def sum_all_int(self, comm_wrapper, value):
        return comm_wrapper.GetWithDeprecationWarning("sum_all_int").SumAll(value)

    def scatter_double(self, comm_wrapper, local_values, source_rank):
        output_list = comm_wrapper.GetWithDeprecationWarning("scatter_double").ScatterDoubles(local_values, source_rank)
        return output_list[0]

    def scatter_int(self, comm_wrapper, local_values, source_rank):
        output_list = comm_wrapper.GetWithDeprecationWarning("scatter_int").ScatterInts(local_values, source_rank)
        return output_list[0]

    def scatterv_double(self, comm_wrapper, local_values, source_rank):
        return comm_wrapper.GetWithDeprecationWarning("scatterv_double").ScattervDoubles(local_values, source_rank)

    def scatterv_int(self, comm_wrapper, local_values, source_rank):
        return comm_wrapper.GetWithDeprecationWarning("scatterv_int").ScattervInts(local_values, source_rank)

    def gather_double(self, comm_wrapper, local_value, gather_rank):
        wraper_list = [local_value,]
        return comm_wrapper.GetWithDeprecationWarning("gather_double").GatherDoubles(wraper_list, gather_rank)

    def gather_int(self, comm_wrapper, local_value, gather_rank):
        wraper_list = [local_value,]
        return comm_wrapper.GetWithDeprecationWarning("gather_int").GatherInts(wraper_list, gather_rank)

    def gatherv_double(self, comm_wrapper, local_values, gather_rank):
        output = comm_wrapper.GetWithDeprecationWarning("gatherv_double").GathervDoubles(local_values, gather_rank)
        if self.__world.Rank() == gather_rank:
            return output
        else:
            return []

    def gatherv_int(self, comm_wrapper, local_values, gather_rank):
        output = comm_wrapper.GetWithDeprecationWarning("gatherv_int").GathervInts(local_values, gather_rank)
        if self.__world.Rank() == gather_rank:
            return output
        else:
            return []

    def allgather_double(self, comm_wrapper, local_value):
        wraper_list = [local_value,]
        return comm_wrapper.GetWithDeprecationWarning("allgather_double").AllGatherDoubles(wraper_list)

    def allgather_int(self, comm_wrapper, local_value):
        wraper_list = [local_value,]
        return comm_wrapper.GetWithDeprecationWarning("allgather_int").AllGatherInts(wraper_list)


