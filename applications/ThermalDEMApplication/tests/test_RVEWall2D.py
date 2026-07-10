import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import auxiliary_functions as Aux

test_folder_name   = 'rve_wall_2d'
solution_file_name = 'reference_solutions.json'

class TestCases(KratosUnittest.TestCase):
    def RVEWall2D_Compression(self):
        parameters_file_name = 'ProjectParametersDEM_Compression.json'
        solution_name        = 'Solution_Compression'
        Aux.CreateAndRunStage(KratosMultiphysics.Model(), test_folder_name, parameters_file_name, solution_file_name, solution_name, 1)

    def RVEWall2D_Expansion(self):
        parameters_file_name = 'ProjectParametersDEM_Expansion.json'
        solution_name        = 'Solution_Expansion'
        Aux.CreateAndRunStage(KratosMultiphysics.Model(), test_folder_name, parameters_file_name, solution_file_name, solution_name, 1)

if __name__ == '__main__':
    KratosUnittest.main()
