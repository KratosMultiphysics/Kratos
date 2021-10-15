# Importing the Kratos Library
import KratosMultiphysics as KM

def GetRankZeroDataCommunicator():
    """In distributed runs, the "Serial" data communicator is defined on all ranks, hence it cannot be used, and we need to create another one that is only defined on rank zero to be used in solvers that do not support MPI
    """
    if KM.IsDistributedRun():
        data_comm_name = "co_simulation_data_comm_rank_zero"
        if not KM.ParallelEnvironment.HasDataCommunicator(data_comm_name):
            from KratosMultiphysics.mpi import DataCommunicatorFactory
            world_data_comm = KM.ParallelEnvironment.GetDataCommunicator("World")
            ranks = [0]
            DataCommunicatorFactory.CreateFromRanksAndRegister(world_data_comm, ranks, data_comm_name)
    else:
        data_comm_name = "Serial"

    return KM.ParallelEnvironment.GetDataCommunicator(data_comm_name)
