try:
    import KratosMultiphysics
    from KratosMultiphysics.mpi import mpi
except ImportError:
    raise Exception("Running in MPI currently requires Kratos-MPI!")

from co_simulation_tools import CoSimulationSpace

class CoSimulationMPISpace(CoSimulationSpace):
    """This requires some MPI-commands exposed to Python
    This is currently only available with Kratos,
    i.e. MPI can only be used with Kratos compiled with MPI
    """

    mpi_op_types_dict = {
        "max" : mpi.MPI_op.MAX,
        "min" : mpi.MPI_op.MIN,
        "sum" : mpi.MPI_op.SUM
    }

    def __init__(self):
        # Precompute rank and size such that they don't have to be recomputed all the time
        self.comm_rank = mpi.rank
        self.comm_size = mpi.size

        if self.comm_size < 2:
            raise Exception("Running in MPI requires at least 2 processes!")

    def IsDistributed(self):
        return True

    def Barrier(self):
        mpi.world.barrier()

    def Rank(self):
        return self.comm_rank

    def Size(self):
        return self.comm_size

    def BroadCast(self, val_to_broadcast, rank_to_broadcast_from):
        self._CheckInputFor_BroadCast(val_to_broadcast, rank_to_broadcast_from)
        return mpi.broadcast(mpi.world, val_to_broadcast, rank_to_broadcast_from)

    def Gather(self, my_value, rank_to_gather_on):
        self._CheckInputFor_Gather(my_value, rank_to_gather_on)
        return mpi.gather(mpi.world, my_value, rank_to_gather_on)

    def AllGather(self, my_value):
        self._CheckInputFor_AllGather(my_value)
        return mpi.allgather(mpi.world, my_value, rank_to_gather_on)

    def GatherV(self, my_values, rank_to_gather_on):
        self._CheckInputFor_GatherV(my_values, rank_to_gather_on)
        return mpi.gatherv(mpi.world, my_values, rank_to_gather_on)

    def Scatter(self, vals_to_scatter, rank_to_scatter_from):
        self._CheckInputFor_Scatter(vals_to_scatter, rank_to_scatter_from)
        # TODO also make use of the int-version
        return mpi.scatter(mpi.world, vals_to_scatter, rank_to_scatter_from)

    def ScatterV(self, vals_to_scatter, rank_to_scatter_from):
        self._CheckInputFor_ScatterV(vals_to_scatter, rank_to_scatter_from)
        # TODO also make use of the int-version
        return mpi.scatterv(mpi.world, vals_to_scatter, rank_to_scatter_from)

    def Reduce(self, my_value, rank_to_reduce_on, mpi_op):
        self._CheckInputFor_Reduce(my_value, rank_to_reduce_on, mpi_op)
        return mpi.reduce(mpi.world, my_value, rank_to_reduce_on, mpi_op_types_dict[mpi_op])

    def AllReduce(self, my_value, rank_to_reduce_on, mpi_op):
        self._CheckInputFor_AllReduce(my_value, rank_to_reduce_on, mpi_op)
        return mpi.allreduce(mpi.world, my_value, mpi_op_types_dict[mpi_op])

