# Importing the Kratos Library
import KratosMultiphysics as KM
if KM.IsDistributedRun():
    from KratosMultiphysics.mpi import DataCommunicatorFactory

def GetRankZeroDataCommunicator():
    """In distributed runs, the "Serial" data communicator is defined on all ranks, hence it cannot be used, and we need to create another one that is only defined on rank zero to be used in solvers that do not support MPI
    """
    if KM.IsDistributedRun():
        data_comm_name = "co_simulation_data_comm_rank_zero"
        if not KM.ParallelEnvironment.HasDataCommunicator(data_comm_name):
            world_data_comm = KM.ParallelEnvironment.GetDataCommunicator("World")
            ranks = [0]
            DataCommunicatorFactory.CreateFromRanksAndRegister(world_data_comm, ranks, data_comm_name)
    else:
        data_comm_name = "Serial"

    return KM.ParallelEnvironment.GetDataCommunicator(data_comm_name)

def CreateDataCommunicatorWithNProcesses(settings):
    if not KM.IsDistributedRun():
        raise Exception("This function can only be used when running with MPI!")

    defaults = KM.Parameters("""{
        "num_processes" : 0,
        "name" : ""
    }""")

    settings.ValidateAndAssignDefaults(defaults)

    num_processes = settings["num_processes"].GetInt()
    name = settings["name"].GetString()

    if num_processes < 1:
        raise Exception('Input for "num_processes" must be > 0, and not {}!'.format(num_processes))
    if name == "":
        raise Exception('Input for "name" cannot be empty!')

    world_data_comm = KM.ParallelEnvironment.GetDataCommunicator("World")
    ranks = list(range(num_processes)) # [ 0 ... num_processes-1 ]
    return DataCommunicatorFactory.CreateFromRanksAndRegister(world_data_comm, ranks, name)
