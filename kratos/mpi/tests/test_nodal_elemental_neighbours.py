import os
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics
from KratosMultiphysics.testing.utilities import ReadModelPart

def GetFilePath(fileName):
    return os.path.dirname(os.path.realpath(__file__)) + "/" + fileName

class TestNodalElementalNeighbours(KratosUnittest.TestCase):

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
        input_filename = GetFilePath("test_mpi_communicator")
        ReadModelPart(input_filename, main_model_part)

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
