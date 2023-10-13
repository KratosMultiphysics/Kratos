# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI

# Other imports
from pathlib import Path
from time import time

class DistributedImportModelPartUtility:

    def __init__(self, main_model_part, settings):
        self.main_model_part = main_model_part
        self.settings = settings
        if settings["model_import_settings"].Has("data_communicator_name"):
            data_comm_name = settings["model_import_settings"]["data_communicator_name"].GetString()
            self.comm = KratosMultiphysics.ParallelEnvironment.GetDataCommunicator(data_comm_name)
        else:
            self.comm = KratosMultiphysics.ParallelEnvironment.GetDefaultDataCommunicator()

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
                "partition_in_memory"                        : false,
                "sub_model_part_list"                        : []
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
                domain_size = 0 # this is not used by Metis, it will be refactored in the future
                verbosity = self.settings["echo_level"].GetInt()

                # Make sure that the condition goes to the same partition as the element is a face of
                sync_conditions = True

                # Original .mdpa file reading
                model_part_io = KratosMultiphysics.ReorderConsecutiveModelPartIO(input_filename, import_flags)

                KratosMultiphysics.Logger.PrintInfo("::[DistributedImportModelPartUtility]::", 'Starting with Metis partitioning of "{}".'.format(input_filename))

                start_time = time()

                if not partition_in_memory:
                    ## Serial partition of the original .mdpa file
                    if self.comm.Rank() == 0:
                        if model_part_import_settings["sub_model_part_list"].size() > 0:
                            no_reorder_model_part_io = KratosMultiphysics.ModelPartIO(input_filename, import_flags)
                            partitioner = KratosMetis.MetisDivideSubModelPartsHeterogeneousInputProcess(no_reorder_model_part_io, model_part_import_settings, number_of_partitions , domain_size, verbosity, sync_conditions)
                        else:
                            partitioner = KratosMetis.MetisDivideHeterogeneousInputProcess(model_part_io, number_of_partitions , domain_size, verbosity, sync_conditions)
                        partitioner.Execute()
                else:
                    # Create a second io that does not reorder the parts while reading from memory
                    serial_model_part_io = KratosMultiphysics.ModelPartIO(input_filename, import_flags)

                    partitioner = KratosMetis.MetisDivideHeterogeneousInputInMemoryProcess(model_part_io, serial_model_part_io, self.comm , domain_size, verbosity, sync_conditions)
                    partitioner.Execute()
                    serial_model_part_io.ReadModelPart(self.main_model_part)

                KratosMultiphysics.Logger.PrintInfo("::[DistributedImportModelPartUtility]::", "Metis divide finished in {0:.{1}f} [s]".format(time()-start_time,1),data_communicator=self.comm)

            else:
                KratosMultiphysics.Logger.PrintInfo("::[DistributedImportModelPartUtility]::", "Metis partitioning not executed.",data_communicator=self.comm)

            self.comm.Barrier()

            ## Reset as input file name the obtained Metis partition one
            if is_single_process_run:
                mpi_input_filename = input_filename
            else:
                base_path = Path(input_filename)
                raw_file_name = base_path.stem
                folder_name = base_path.parent / Path(str(raw_file_name) + "_partitioned")

                mpi_input_filename = str(folder_name / Path(str(raw_file_name) + "_"+str(self.comm.Rank())))

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
        ParallelFillCommunicator = KratosMPI.ParallelFillCommunicator(self.main_model_part.GetRootModelPart(), self.comm)
        ParallelFillCommunicator.Execute()

        KratosMultiphysics.Logger.PrintInfo("::[DistributedImportModelPartUtility]::", "MPI communicators constructed.")
