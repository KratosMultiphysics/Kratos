from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI

import os

class DistributedImportModelPartUtility(object):

    def __init__(self, main_model_part, settings):
        self.main_model_part = main_model_part
        self.settings = settings
        self.comm = KratosMultiphysics.DataCommunicator.GetDefault()

    def ExecutePartitioningAndReading(self):
        warning_msg  = 'Calling "ExecutePartitioningAndReading" which is DEPRECATED\n'
        warning_msg += 'Please use "ImportModelPart" instead'
        KratosMultiphysics.Logger.PrintWarning("::[DistributedImportModelPartUtility]::", warning_msg, data_communicator=self.comm)

        self.ImportModelPart()

    def ImportModelPart(self):
        model_part_import_settings = self.settings["model_import_settings"]
        input_type = model_part_import_settings["input_type"].GetString()

        # in single process runs, do not call metis (no partitioning is necessary)
        is_single_process_run = (self.comm.Size() == 1)

        if input_type == "mdpa":
            default_settings = KratosMultiphysics.Parameters("""{
                "input_filename"                             : "",
                "skip_timer"                                 : true,
                "ignore_variables_not_in_solution_step_data" : false,
                "perform_partitioning"                       : true,
                "partition_in_memory"                        : false
            }""")

            # cannot validate as this might contain other settings too
            model_part_import_settings.AddMissingParameters(default_settings)

            input_filename = model_part_import_settings["input_filename"].GetString()

            perform_partitioning = model_part_import_settings["perform_partitioning"].GetBool()
            partition_in_memory = model_part_import_settings["partition_in_memory"].GetBool()

            # Setting some mdpa-import-related flags
            import_flags = KratosMultiphysics.ModelPartIO.READ

            if model_part_import_settings["skip_timer"].GetBool():
                import_flags = KratosMultiphysics.ModelPartIO.SKIP_TIMER|import_flags

            if model_part_import_settings["ignore_variables_not_in_solution_step_data"].GetBool():
                import_flags = KratosMultiphysics.ModelPartIO.IGNORE_VARIABLES_ERROR|import_flags

            if not is_single_process_run and perform_partitioning:
                import KratosMultiphysics.MetisApplication as KratosMetis

                # Partition of the original .mdpa file
                number_of_partitions = self.comm.Size() # Number of partitions equals the number of processors
                domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
                verbosity = self.settings["echo_level"].GetInt()

                # Make sure that the condition goes to the same partition as the element is a face of
                sync_conditions = True

                # Original .mdpa file reading
                model_part_io = KratosMultiphysics.ReorderConsecutiveModelPartIO(input_filename)

                if not partition_in_memory:
                    ## Serial partition of the original .mdpa file
                    if self.comm.Rank() == 0:
                        partitioner = KratosMetis.MetisDivideHeterogeneousInputProcess(model_part_io, number_of_partitions , domain_size, verbosity, sync_conditions)
                        partitioner.Execute()

                        KratosMultiphysics.Logger.PrintInfo("::[DistributedImportModelPartUtility]::", "Metis divide finished.")
                else:
                    # Create a second io that does not reorder the parts while reading from memory
                    serial_model_part_io = KratosMultiphysics.ModelPartIO(input_filename, import_flags)

                    partitioner = KratosMetis.MetisDivideHeterogeneousInputInMemoryProcess(model_part_io, serial_model_part_io, number_of_partitions , domain_size, verbosity, sync_conditions)
                    partitioner.Execute()
                    serial_model_part_io.ReadModelPart(self.main_model_part)

                    KratosMultiphysics.Logger.PrintInfo("::[DistributedImportModelPartUtility]::", "Metis divide finished.",data_communicator=self.comm)

            else:
                KratosMultiphysics.Logger.PrintInfo("::[DistributedImportModelPartUtility]::", "Metis partitioning not executed.",data_communicator=self.comm)

            self.comm.Barrier()

            ## Reset as input file name the obtained Metis partition one
            if is_single_process_run:
                mpi_input_filename = input_filename
            else:
                mpi_input_filename = os.path.join(input_filename+"_partitioned", input_filename + "_" + str(self.comm.Rank())
            model_part_import_settings["input_filename"].SetString(mpi_input_filename)

            ## Read the new generated *.mdpa files
            if not partition_in_memory or is_single_process_run:
                KratosMultiphysics.ModelPartIO(mpi_input_filename, import_flags).ReadModelPart(self.main_model_part)

        elif input_type == "rest":
            from KratosMultiphysics.mpi.distributed_restart_utility import DistributedRestartUtility as RestartUtility
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
        ParallelFillCommunicator = KratosMPI.ParallelFillCommunicator(self.main_model_part.GetRootModelPart())
        ParallelFillCommunicator.Execute()

        KratosMultiphysics.Logger.PrintInfo("::[DistributedImportModelPartUtility]::", "MPI communicators constructed.")
