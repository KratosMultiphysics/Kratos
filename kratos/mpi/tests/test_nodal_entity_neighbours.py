import os
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics
from KratosMultiphysics.testing.utilities import ReadModelPart

def GetFilePath(fileName):
    return os.path.dirname(os.path.realpath(__file__)) + "/" + fileName

class TestNodalEntityNeighbours(KratosUnittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.current_model = KratosMultiphysics.Model()
        cls.main_model_part = cls.current_model.CreateModelPart("main_model_part")
        cls.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]=2

        ## Add variables to the model part
        cls.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        cls.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        cls.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        cls.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)

        ## Serial partition of the original .mdpa file
        input_filename = GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_mpi_communicator")
        ReadModelPart(input_filename, cls.main_model_part)

    def test_global_elemental_neighbour_pointers(self):
        #compute nodal neighbours
        neighbour_finder = KratosMultiphysics.FindGlobalNodalElementalNeighboursProcess(self.main_model_part)
        neighbour_finder.Execute()

        #obtain the ids of all the neighbours
        found_ids = neighbour_finder.GetNeighbourIds(self.main_model_part.Nodes)

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

        self.__CheckNeighbourIds(found_ids, reference_ids)

    def test_global_condition_neighbour_pointers(self):
        #compute nodal neighbours
        neighbour_finder = KratosMultiphysics.FindGlobalNodalConditionNeighboursProcess(self.main_model_part)
        neighbour_finder.Execute()

        #obtain the ids of all the neighbours
        found_ids = neighbour_finder.GetNeighbourIds(self.main_model_part.Nodes)

        reference_ids = {
            9: [4, 5],
            8: [5, 6],
            7: [6, 7],
            6: [3, 4],
            5: [],
            4: [7, 8],
            3: [2, 3],
            2: [1, 2],
            1: [1, 8]
        }

        self.__CheckNeighbourIds(found_ids, reference_ids)

    def __CheckNeighbourIds(self, found_ids, reference_ids):
        for key,values in found_ids.items():
            ref_values = reference_ids[key]
            self.assertEqual(len(ref_values), len(values))
            for i in range(len(ref_values)):
                self.assertEqual(values[i],ref_values[i])

if __name__ == '__main__':
    KratosUnittest.main()
