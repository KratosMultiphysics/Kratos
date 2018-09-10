from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("TrilinosApplication")

# Import applications
import KratosMultiphysics.TrilinosApplication as KratosTrilinos

# Other imports
import restart_utility

class TrilinosRestartUtility(restart_utility.RestartUtility):
    """
    This class overwrites the methods that are different
    in MPI parallel execution

    See restart_utility.py for more information.
    """
    def __init__(self, model_part, settings):
        # Construct the base class
        super(TrilinosRestartUtility, self).__init__(model_part, settings)

        self.set_mpi_communicator = settings["set_mpi_communicator"].GetBool()

    #### Protected functions ####

    def _GetFileLabelLoad(self):
        return str(KratosMPI.mpi.rank) + '_' + self.input_file_label

    def _GetFileLabelSave(self, file_label):
        return str(KratosMPI.mpi.rank) + '_' + str(file_label)

    def _ExecuteAfterLoad(self):
        if self.set_mpi_communicator:
            KratosTrilinos.ParallelFillCommunicator(self.main_model_part.GetRootModelPart()).Execute()

    def _PrintOnRankZero(self, *args):
        KratosMPI.mpi.world.barrier()
        if KratosMPI.mpi.rank == 0:
            KratosMultiphysics.Logger.PrintInfo(" ".join(map(str,args)))
