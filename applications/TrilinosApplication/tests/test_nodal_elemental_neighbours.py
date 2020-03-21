from __future__ import print_function, absolute_import, division

import os
import sys
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics
import KratosMultiphysics.kratos_utilities as kratos_utilities
from KratosMultiphysics.mpi import distributed_import_model_part_utility

def GetFilePath(fileName):
    return os.path.dirname(os.path.realpath(__file__)) + "/" + fileName

class TestNodalElementalNeighbours(KratosUnittest.TestCase):

    def tearDown(self):
        kratos_comm  = KratosMultiphysics.DataCommunicator.GetDefault()
        rank = kratos_comm.Rank()
        kratos_utilities.DeleteFileIfExisting("test_mpi_communicator.time")
        kratos_utilities.DeleteFileIfExisting("test_mpi_communicator_"+str(rank)+".mdpa")
        kratos_utilities.DeleteFileIfExisting("test_mpi_communicator_"+str(rank)+".time")
        kratos_comm.Barrier()

    def _ReadSerialModelPart(self,model_part, mdpa_file_name):
        import_flags = KratosMultiphysics.ModelPartIO.READ | KratosMultiphysics.ModelPartIO.SKIP_TIMER
        KratosMultiphysics.ModelPartIO(mdpa_file_name, import_flags).ReadModelPart(model_part)

    def _ReadDistributedModelPart(self,model_part, mdpa_file_name):

        importer_settings = KratosMultiphysics.Parameters("""{
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": \"""" + mdpa_file_name + """\",
                "partition_in_memory" : false
            },
            "echo_level" : 0
        }""")

        model_part_import_util = distributed_import_model_part_utility.DistributedImportModelPartUtility(model_part, importer_settings)
        model_part_import_util.ImportModelPart()
        model_part_import_util.CreateCommunicators()

    def _ReadModelPart(self, input_filename, mp):
        kratos_comm  = KratosMultiphysics.DataCommunicator.GetDefault()

        if kratos_comm.IsDistributed():
            self._ReadDistributedModelPart(mp, input_filename)
        else:
            self._ReadSerialModelPart(mp, input_filename)

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
        neighbour_finder = KratosMultiphysics.FindGlobalNodalElementalNeighboursProcess(kratos_comm, main_model_part)
        neighbour_finder.Execute()

        #obtain the ids of all the neighbours
        found_ids = neighbour_finder.GetNeighbourIds(main_model_part.Nodes)

        reference_ids = {
            1 : [1,2],
            2 : [2,3,4],
            3 : [4],
            4 : [1,5,6],
            5 : [1,2,3,6,7,8],
            6 : [3,4,8],
            7 : [5],
            8 : [5,6,7],
            9 : [7,8]
        }

        #do the check
        for key,values in found_ids.items():
            ref_values = reference_ids[key]
            for i in range(len(ref_values)):
                self.assertEqual(values[i],ref_values[i])

if __name__ == '__main__':
    KratosUnittest.main()
