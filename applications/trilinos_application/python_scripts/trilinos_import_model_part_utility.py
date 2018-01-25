from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("MetisApplication","TrilinosApplication")

# Import applications
import KratosMultiphysics.MetisApplication as KratosMetis
import KratosMultiphysics.TrilinosApplication as KratosTrilinos


class TrilinosImportModelPartUtility():

    def __init__(self, main_model_part, settings):
        self.main_model_part = main_model_part
        self.settings = settings

        # Otherwise stated, always perform the Metis partitioning
        if (self.settings["model_import_settings"].Has("perform_partitioning") == False):
            self.settings["model_import_settings"].AddEmptyValue("perform_partitioning")
            self.settings["model_import_settings"]["perform_partitioning"].SetBool(True)

        if KratosMPI.mpi.size < 2:
            raise NameError("MPI number of processors is 1.")

    def ExecutePartitioningAndReading(self):
        if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):
            input_filename = self.settings["model_import_settings"]["input_filename"].GetString()
            perform_partitioning = self.settings["model_import_settings"]["perform_partitioning"].GetBool()
            ## Serial partition of the original .mdpa file
            if ((KratosMPI.mpi.rank == 0) and (perform_partitioning == True)):

                # Original .mdpa file reading
                model_part_io = KratosMultiphysics.ReorderConsecutiveModelPartIO(input_filename)

                # Partition of the original .mdpa file
                number_of_partitions = KratosMPI.mpi.size # Number of partitions equals the number of processors
                domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
                verbosity = 1
                sync_conditions = True # Make sure that the condition goes to the same partition as the element is a face of
                partitioner = KratosMetis.MetisDivideHeterogeneousInputProcess(model_part_io, number_of_partitions , domain_size, verbosity, sync_conditions)
                partitioner.Execute()

                print("Metis divide finished.")

            else:
                if (KratosMPI.mpi.rank == 0):
                    print("Metis partitioning not executed.")

            KratosMPI.mpi.world.barrier()

            ## Reset as input file name the obtained Metis partition one
            mpi_input_filename = input_filename + "_" + str(KratosMPI.mpi.rank)
            self.settings["model_import_settings"]["input_filename"].SetString(mpi_input_filename)

            ## Read the new generated *.mdpa files
            KratosMultiphysics.ModelPartIO(mpi_input_filename).ReadModelPart(self.main_model_part)

        else:
            raise Exception("Other input options are not yet implemented.")

    def CreateCommunicators(self):
        ## Construct and execute the MPICommunicator
        KratosMetis.SetMPICommunicatorProcess(self.main_model_part).Execute()

        ## Construct and execute the Parallel fill communicator
        ParallelFillCommunicator = KratosTrilinos.ParallelFillCommunicator(self.main_model_part.GetRootModelPart())
        ParallelFillCommunicator.Execute()

        if KratosMPI.mpi.rank == 0 :
            print("MPI communicators constructed.")
