from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI
from compare_two_files_check_process import CompareTwoFilesCheckProcess

def Factory(settings, current_model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return CompareTwoFilesCheckProcessMPI(settings["Parameters"])

class CompareTwoFilesCheckProcessMPI(CompareTwoFilesCheckProcess):
    """Inserts MPI barrier before check to ensure results files are present."""

    def ExecuteFinalize(self):
        KratosMPI.mpi.world.barrier()
        super(CompareTwoFilesCheckProcessMPI, self).ExecuteFinalize()
