import os
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics
from KratosMultiphysics.testing.utilities import ReadModelPart

def GetFilePath(fileName):
    return os.path.dirname(os.path.realpath(__file__)) + "/" + fileName

class TestNeighbours(KratosUnittest.TestCase):

    def test_global_neighbour_pointers(self):

        kratos_comm  = KratosMultiphysics.DataCommunicator.GetDefault()

        current_model = KratosMultiphysics.Model()
        main_model_part = current_model.CreateModelPart("main_model_part")
        main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]=2

        ## Add variables to the model part
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)

        ## Serial partition of the original .mdpa file
        input_filename = GetFilePath("test_files/mdpa_files/test_mpi_communicator")
        ReadModelPart(input_filename, main_model_part)

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

        for key, values in found_ids.items():
            ref_values = reference_ids[key]
            self.assertListEqual(values, ref_values)


if __name__ == '__main__':
    KratosUnittest.main()
