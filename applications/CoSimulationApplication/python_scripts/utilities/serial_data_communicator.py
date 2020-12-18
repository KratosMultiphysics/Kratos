# Importing the Kratos Library
import KratosMultiphysics as KM
default_data_comm = KM.DataCommunicator.GetDefault()

class SerialDataCommunicator:
    """Auxiliar object for MPI runs in which also serial solvers participate
    The serial data communicator cannot be used as it also exists on every rank
    (see parallel_environment.cpp)
    """
    def __init__(self):
        self.rank = default_data_comm.Rank()

    def IsDefinedOnThisRank(self):
        return self.rank == 0
