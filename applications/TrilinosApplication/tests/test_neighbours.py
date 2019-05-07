from __future__ import print_function, absolute_import, division

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI
import KratosMultiphysics.MetisApplication as KratosMetis
import KratosMultiphysics.TrilinosApplication as TrilinosApplication
import KratosMultiphysics.kratos_utilities as kratos_utilities

import trilinos_import_model_part_utility
import os



def GetFilePath(fileName):
    return os.path.dirname(os.path.realpath(__file__)) + "/" + fileName


class TestMPICommunicator(KratosUnittest.TestCase):

    def tearDown(self):
        kratos_comm  = KratosMultiphysics.DataCommunicator.GetDefault()
        rank = kratos_comm.Rank()
        if rank == 0:
            kratos_utilities.DeleteFileIfExisting("test_mpi_communicator.time")
        kratos_utilities.DeleteFileIfExisting("test_mpi_communicator_"+str(rank)+".mdpa")
        kratos_utilities.DeleteFileIfExisting("test_mpi_communicator_"+str(rank)+".time")
        kratos_comm.Barrier()

    def _ReadModelPart(self, input_filename, mp):
        kratos_comm  = KratosMultiphysics.DataCommunicator.GetDefault()

        if(kratos_comm.IsDistributed()):
            # import_settings = KratosMultiphysics.Parameters("""{
            #     "echo_level" : 0,
            #     "model_import_settings" : {
            #         "input_type" : "mdpa",
            #         "input_filename" :""
            #     }
            # } """)
            # import_settings["model_import_settings"]["input_filename"].SetString(input_filename)

            # TrilinosModelPartImporter = trilinos_import_model_part_utility.TrilinosImportModelPartUtility(mp, import_settings)
            # TrilinosModelPartImporter.ImportModelPart()
            # TrilinosModelPartImporter.CreateCommunicators()
            # print(mp)
            if kratos_comm.Rank() == 0 :

                # Original .mdpa file reading
                model_part_io = KratosMultiphysics.ModelPartIO(input_filename)

                # Partition of the original .mdpa file
                number_of_partitions = kratos_comm.Size() # Number of partitions equals the number of processors
                domain_size = mp.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
                verbosity = 0
                sync_conditions = True # Make sure that the condition goes to the same partition as the element is a face of

                partitioner = KratosMetis.MetisDivideHeterogeneousInputProcess(model_part_io, number_of_partitions , domain_size, verbosity, sync_conditions)
                partitioner.Execute()

            kratos_comm.Barrier()

            # Read the partitioned .mdpa files
            mpi_input_filename = input_filename + "_" + str(kratos_comm.Rank())
            model_part_io = KratosMultiphysics.ModelPartIO(mpi_input_filename)
            model_part_io.ReadModelPart(mp)
            print(mp)
        else:
            model_part_io = KratosMultiphysics.ModelPartIO(input_filename)
            model_part_io.ReadModelPart(mp)

        print("reading finished")

    def test_global_neighbour_pointers(self):

        kratos_comm  = KratosMultiphysics.DataCommunicator.GetDefault()

        current_model = KratosMultiphysics.Model()
        main_model_part = current_model.CreateModelPart("main_model_part")

        ## Add variables to the model part
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)

        ## Serial partition of the original .mdpa file
        input_filename = GetFilePath("test_mpi_communicator")
        self._ReadModelPart(input_filename, main_model_part)

        #compute nodal neighbours
        neighbour_finder = KratosMultiphysics.FindGlobalNodalNeighboursProcess(kratos_comm, main_model_part)
        neighbour_finder.Execute()

        #obtain the ids of all the neighbours
        found_ids = neighbour_finder.GetNeighbourIds(kratos_comm, main_model_part.Nodes)


        reference_ids = {
            1 : [2,4,5],
            2 : [1,3,5,6],
            3 : [2,6],
            4 : [1,5,7,8],
            5 : [1,2,4,6,8,9],
            6 : [2,3,5,9],
            7 : [4,8],
            8 : [4,5,7,9],
            9 : [5,6,8]
        }

        #do the check
        for key,values in found_ids.items():
            ref_values = reference_ids[key]
            for i in range(len(values)):
                self.assertEqual(values[i],ref_values[i])

if __name__ == '__main__':
    KratosUnittest.main()
