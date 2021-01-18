import os
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics
from KratosMultiphysics.kratos_utilities import DeleteDirectoryIfExisting
from KratosMultiphysics.mpi import distributed_import_model_part_utility

def GetFilePath(fileName):
    return os.path.dirname(os.path.realpath(__file__)) + "/" + fileName

class TestNeighbours(KratosUnittest.TestCase):

    def tearDown(self):
        DeleteDirectoryIfExisting("test_mpi_communicator_partitioned")

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
        neighbour_finder = KratosMultiphysics.FindGlobalNodalNeighboursProcess(kratos_comm, main_model_part)
        neighbour_finder.Execute()

        #obtain the ids of all the neighbours
        found_ids = neighbour_finder.GetNeighbourIds(main_model_part.Nodes)

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
