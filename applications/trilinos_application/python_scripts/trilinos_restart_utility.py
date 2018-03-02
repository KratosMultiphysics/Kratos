from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("MetisApplication","TrilinosApplication")

# Import applications
import KratosMultiphysics.MetisApplication as KratosMetis
import KratosMultiphysics.TrilinosApplication as KratosTrilinos

# Other imports
import restart_utility
import os

class TrilinosRestartUtility(RestartUtility):
    def __init__(self, model_part, settings):
        # Construct the base class
        super(TrilinosRestartUtility, self).__init__(model_part, settings)

    def __GetFileNameLoad(self):
        problem_path = os.getcwd()
        return os.path.join(problem_path, self.input_filename + '_' + str(KratosMPI.mpi.rank) + "_" + str(self.input_file_label))

    def __GetFileNameSave(self, time):
        return self.input_filename + '_' + str(KratosMPI.mpi.rank) + '_' + str(time)

    def __ExecuteAfterLaod(self):
        import trilinos_import_model_part_utility
        trilinos_model_part_importer = trilinos_import_model_part_utility.TrilinosImportModelPartUtility(self.model_part, Parameters('''{"model_import_settings":[]}'''))
        trilinos_model_part_importer.CreateCommunicators() # parallel fill communicator

    def __PrintOnRankZero(self, *args):
        KratosMPI.mpi.world.barrier()
        if KratosMPI.mpi.rank == 0:
            KratosMultiphysics.Logger.PrintInfo(" ".join(map(str,args)))