from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("MetisApplication","TrilinosApplication")

# Import applications
import KratosMultiphysics.TrilinosApplication as KratosTrilinos


class TrilinosImportModelPartUtility():

    def __init__(self, main_model_part, settings):
        self.main_model_part = main_model_part
        self.settings = settings

        if KratosMPI.mpi.size < 2:
            raise NameError("MPI number of processors is 1.")

    def ExecutePartitioningAndReading(self):
        warning_msg  = 'Calling "ExecutePartitioningAndReading" which is DEPRECATED\n'
        warning_msg += 'Please use "ImportModelPart" instead'
        if (KratosMPI.mpi.rank == 0):
            KratosMultiphysics.Logger.PrintWarning("TrilinosImportModelPartUtility", warning_msg)

        self.ImportModelPart()

    def ImportModelPart(self):
        input_type = self.settings["model_import_settings"]["input_type"].GetString()

        if input_type == "mdpa":
            input_filename = self.settings["model_import_settings"]["input_filename"].GetString()

            # Unless otherwise stated, always perform the Metis partitioning
            if not self.settings["model_import_settings"].Has("perform_partitioning"):
                self.settings["model_import_settings"].AddEmptyValue("perform_partitioning")
                self.settings["model_import_settings"]["perform_partitioning"].SetBool(True)

            perform_partitioning = self.settings["model_import_settings"]["perform_partitioning"].GetBool()

            if perform_partitioning == True:
                KratosMultiphysics.CheckRegisteredApplications("MetisApplication")
                import KratosMultiphysics.MetisApplication as KratosMetis

                ## Serial partition of the original .mdpa file
                if KratosMPI.mpi.rank == 0:

                    # Original .mdpa file reading
                    model_part_io = KratosMultiphysics.ReorderConsecutiveModelPartIO(input_filename)

                    # Partition of the original .mdpa file
                    number_of_partitions = KratosMPI.mpi.size # Number of partitions equals the number of processors
                    domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
                    verbosity = self.settings["echo_level"].GetInt()
                    sync_conditions = True # Make sure that the condition goes to the same partition as the element is a face of
                    partitioner = KratosMetis.MetisDivideHeterogeneousInputProcess(model_part_io, number_of_partitions , domain_size, verbosity, sync_conditions)
                    partitioner.Execute()

                    KratosMultiphysics.Logger.PrintInfo("::[TrilinosImportModelPartUtility]::", "Metis divide finished.")

            else:
                if (KratosMPI.mpi.rank == 0):
                    KratosMultiphysics.Logger.PrintInfo("::[TrilinosImportModelPartUtility]::", "Metis partitioning not executed.")

            KratosMPI.mpi.world.barrier()

            ## Reset as input file name the obtained Metis partition one
            mpi_input_filename = input_filename + "_" + str(KratosMPI.mpi.rank)
            self.settings["model_import_settings"]["input_filename"].SetString(mpi_input_filename)

            ## Read the new generated *.mdpa files
            KratosMultiphysics.ModelPartIO(mpi_input_filename).ReadModelPart(self.main_model_part)

        elif input_type == "rest":
            from trilinos_restart_utility import TrilinosRestartUtility as RestartUtility
            restart_settings = self.settings["model_import_settings"].Clone()

            restart_settings.RemoveValue("input_type")

            restart_settings.AddEmptyValue("set_mpi_communicator")
            restart_settings["set_mpi_communicator"].SetBool(False)

            if not restart_settings.Has("restart_load_file_label"):
                raise Exception('"restart_load_file_label" must be specified when starting from a restart-file!')

            RestartUtility(self.main_model_part, restart_settings).LoadRestart()

        elif input_type == "use_input_model_part":
            pass

        else:
            raise Exception("Other input options are not yet implemented.")

    def CreateCommunicators(self):
        ## Construct and execute the Parallel fill communicator (also sets the MPICommunicator)
        ParallelFillCommunicator = KratosTrilinos.ParallelFillCommunicator(self.main_model_part.GetRootModelPart())
        ParallelFillCommunicator.Execute()

        if KratosMPI.mpi.rank == 0 :
            KratosMultiphysics.Logger.PrintInfo("::[TrilinosImportModelPartUtility]::", "MPI communicators constructed.")
