# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosCoSimulationApplication import *
application = KratosCoSimulationApplication()
application_name = "KratosCoSimulationApplication"

_ImportApplication(application, application_name)

def __ModuleInitDetail():
    """
    Create a DataCommunicator that only contains rank zero and is undefined in other ranks
    This is necessary for solvers that can only run in serial
    It is defined as a function to avoid polluting the Kratos namespace with local variables.
    """
    import KratosMultiphysics as KM
    if KM.IsDistributedRun():
        from KratosMultiphysics.mpi import DataCommunicatorFactory
        data_comm_name = "co_simulation_data_comm_rank_zero"
        world_data_comm = KM.ParallelEnvironment.GetDataCommunicator("World")
        ranks = [0]
        DataCommunicatorFactory.CreateFromRanksAndRegister(world_data_comm, ranks, data_comm_name)

__ModuleInitDetail()
