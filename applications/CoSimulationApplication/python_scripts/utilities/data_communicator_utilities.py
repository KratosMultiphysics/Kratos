# Importing the Kratos Library
import KratosMultiphysics as KM
if KM.IsDistributedRun():
    from KratosMultiphysics.mpi import DataCommunicatorFactory

def GetRankZeroDataCommunicator():
    """In distributed runs, the "Serial" data communicator is defined on all ranks, hence it cannot be used, and we need to create another one that is only defined on rank zero to be used in solvers that do not support MPI
    """
    if KM.IsDistributedRun():
        data_comm_name = "co_simulation_data_comm_rank_zero"
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

    world_data_comm = KM.ParallelEnvironment.GetDataCommunicator("World")

    if num_processes < 1:
        raise Exception('Input for "num_processes" ({}) must be > 0!'.format(num_processes))
    if num_processes > world_data_comm.Size():
        raise Exception('Input for "num_processes" ({}) cannot be larger than the number of MPI processes ({})!'.format(num_processes, world_data_comm.Size()))
    if name == "":
        raise Exception('Input for "name" cannot be empty!')

    ranks = list(range(num_processes)) # [ 0 ... num_processes-1 ]
    return DataCommunicatorFactory.CreateFromRanksAndRegister(world_data_comm, ranks, name)
