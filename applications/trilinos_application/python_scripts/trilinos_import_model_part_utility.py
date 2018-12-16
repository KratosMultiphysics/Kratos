from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI

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
        model_part_import_settings = self.settings["model_import_settings"]
        input_type = model_part_import_settings["input_type"].GetString()

        if input_type == "mdpa":
            input_filename = model_part_import_settings["input_filename"].GetString()

            # Unless otherwise stated, always perform the Metis partitioning
            if not model_part_import_settings.Has("perform_partitioning"):
                model_part_import_settings.AddEmptyValue("perform_partitioning")
                model_part_import_settings["perform_partitioning"].SetBool(True)

            perform_partitioning = model_part_import_settings["perform_partitioning"].GetBool()

            # Setting some mdpa-import-related flags
            import_flags = KratosMultiphysics.ModelPartIO.READ
            if model_part_import_settings.Has("ignore_variables_not_in_solution_step_data"):
                if model_part_import_settings["ignore_variables_not_in_solution_step_data"].GetBool():
                    import_flags = KratosMultiphysics.ModelPartIO.IGNORE_VARIABLES_ERROR|import_flags
            skip_timer = True
            if model_part_import_settings.Has("skip_timer"):
                skip_timer = model_part_import_settings["skip_timer"].GetBool()
            if skip_timer:
                import_flags = KratosMultiphysics.ModelPartIO.SKIP_TIMER|import_flags

            # Select the partitioning method (File by default)
            partition_in_memory = False
            if model_part_import_settings.Has("partition_in_memory"):
                partition_in_memory = model_part_import_settings["partition_in_memory"].GetBool()

            if perform_partitioning == True:
                import KratosMultiphysics.MetisApplication as KratosMetis

                # Partition of the original .mdpa file
                number_of_partitions = KratosMPI.mpi.size # Number of partitions equals the number of processors
                domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
                verbosity = self.settings["echo_level"].GetInt()

                # Make sure that the condition goes to the same partition as the element is a face of
                sync_conditions = True

                # Original .mdpa file reading
                model_part_io = KratosMultiphysics.ReorderConsecutiveModelPartIO(input_filename)

                if not partition_in_memory:
                    ## Serial partition of the original .mdpa file
                    if KratosMPI.mpi.rank == 0:
                        partitioner = KratosMetis.MetisDivideHeterogeneousInputProcess(model_part_io, number_of_partitions , domain_size, verbosity, sync_conditions)
                        partitioner.Execute()

                        KratosMultiphysics.Logger.PrintInfo("::[TrilinosImportModelPartUtility]::", "Metis divide finished.")
                else:
                    # Create a second io that does not reorder the parts while reading from memory
                    serial_model_part_io = KratosMultiphysics.ModelPartIO(input_filename, import_flags)

                    partitioner = KratosMetis.MetisDivideHeterogeneousInputInMemoryProcess(model_part_io, serial_model_part_io, number_of_partitions , domain_size, verbosity, sync_conditions)
                    partitioner.Execute()
                    serial_model_part_io.ReadModelPart(self.main_model_part)

                    if KratosMPI.mpi.rank == 0:
                        KratosMultiphysics.Logger.PrintInfo("::[TrilinosImportModelPartUtility]::", "Metis divide finished.")

            else:
                if (KratosMPI.mpi.rank == 0):
                    KratosMultiphysics.Logger.PrintInfo("::[TrilinosImportModelPartUtility]::", "Metis partitioning not executed.")

            KratosMPI.mpi.world.barrier()

            ## Reset as input file name the obtained Metis partition one
            mpi_input_filename = input_filename + "_" + str(KratosMPI.mpi.rank)
            model_part_import_settings["input_filename"].SetString(mpi_input_filename)

            ## Read the new generated *.mdpa files
            if not partition_in_memory:
                KratosMultiphysics.ModelPartIO(mpi_input_filename, import_flags).ReadModelPart(self.main_model_part)

        elif input_type == "rest":
            from trilinos_restart_utility import TrilinosRestartUtility as RestartUtility
            restart_settings = model_part_import_settings.Clone()

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