import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import auxiliary_functions as Aux

test_folder_name   = 'colliding_spheres_3d'
solution_file_name = 'reference_solutions.json'

class TestCases(KratosUnittest.TestCase):
    def CollidingSpheres3D_Mixed_01(self):
        parameters_file_name = 'ProjectParametersDEM_Mixed_01.json'
        solution_name        = 'Solution_Mixed_01'
        Aux.CreateAndRunStage(KratosMultiphysics.Model(), test_folder_name, parameters_file_name, solution_file_name, solution_name, 1)

    def CollidingSpheres3D_ConductionDir_01(self):
        parameters_file_name = 'ProjectParametersDEM_ConductionDir_01.json'
        solution_name        = 'Solution_ConductionDir_01'
        Aux.CreateAndRunStage(KratosMultiphysics.Model(), test_folder_name, parameters_file_name, solution_file_name, solution_name, 1)

    def CollidingSpheres3D_ConductionDir_02(self):
        parameters_file_name = 'ProjectParametersDEM_ConductionDir_02.json'
        solution_name        = 'Solution_ConductionDir_02'
        Aux.CreateAndRunStage(KratosMultiphysics.Model(), test_folder_name, parameters_file_name, solution_file_name, solution_name, 1)

    def CollidingSpheres3D_ConductionDir_03(self):
        parameters_file_name = 'ProjectParametersDEM_ConductionDir_03.json'
        solution_name        = 'Solution_ConductionDir_03'
        Aux.CreateAndRunStage(KratosMultiphysics.Model(), test_folder_name, parameters_file_name, solution_file_name, solution_name, 1)

    def CollidingSpheres3D_ConductionDir_04(self):
        parameters_file_name = 'ProjectParametersDEM_ConductionDir_04.json'
        solution_name        = 'Solution_ConductionDir_04'
        Aux.CreateAndRunStage(KratosMultiphysics.Model(), test_folder_name, parameters_file_name, solution_file_name, solution_name, 1)

    def CollidingSpheres3D_ConductionDir_05(self):
        parameters_file_name = 'ProjectParametersDEM_ConductionDir_05.json'
        solution_name        = 'Solution_ConductionDir_05'
        Aux.CreateAndRunStage(KratosMultiphysics.Model(), test_folder_name, parameters_file_name, solution_file_name, solution_name, 1)

    def CollidingSpheres3D_ConductionIndir_01(self):
        parameters_file_name = 'ProjectParametersDEM_ConductionIndir_01.json'
        solution_name        = 'Solution_ConductionIndir_01'
        Aux.CreateAndRunStage(KratosMultiphysics.Model(), test_folder_name, parameters_file_name, solution_file_name, solution_name, 1)

    def CollidingSpheres3D_ConductionIndir_02(self):
        parameters_file_name = 'ProjectParametersDEM_ConductionIndir_02.json'
        solution_name        = 'Solution_ConductionIndir_02'
        Aux.CreateAndRunStage(KratosMultiphysics.Model(), test_folder_name, parameters_file_name, solution_file_name, solution_name, 1)

    def CollidingSpheres3D_ConductionIndir_03(self):
        parameters_file_name = 'ProjectParametersDEM_ConductionIndir_03.json'
        solution_name        = 'Solution_ConductionIndir_03'
        Aux.CreateAndRunStage(KratosMultiphysics.Model(), test_folder_name, parameters_file_name, solution_file_name, solution_name, 1)

    def CollidingSpheres3D_ConductionIndir_04(self):
        parameters_file_name = 'ProjectParametersDEM_ConductionIndir_04.json'
        solution_name        = 'Solution_ConductionIndir_04'
        Aux.CreateAndRunStage(KratosMultiphysics.Model(), test_folder_name, parameters_file_name, solution_file_name, solution_name, 1)

    def CollidingSpheres3D_Convection_01(self):
        parameters_file_name = 'ProjectParametersDEM_Convection_01.json'
        solution_name        = 'Solution_Convection_01'
        Aux.CreateAndRunStage(KratosMultiphysics.Model(), test_folder_name, parameters_file_name, solution_file_name, solution_name, 1)

    def CollidingSpheres3D_Convection_02(self):
        parameters_file_name = 'ProjectParametersDEM_Convection_02.json'
        solution_name        = 'Solution_Convection_02'
        Aux.CreateAndRunStage(KratosMultiphysics.Model(), test_folder_name, parameters_file_name, solution_file_name, solution_name, 1)

    def CollidingSpheres3D_Convection_03(self):
        parameters_file_name = 'ProjectParametersDEM_Convection_03.json'
        solution_name        = 'Solution_Convection_03'
        Aux.CreateAndRunStage(KratosMultiphysics.Model(), test_folder_name, parameters_file_name, solution_file_name, solution_name, 1)

    def CollidingSpheres3D_Convection_04(self):
        parameters_file_name = 'ProjectParametersDEM_Convection_04.json'
        solution_name        = 'Solution_Convection_04'
        Aux.CreateAndRunStage(KratosMultiphysics.Model(), test_folder_name, parameters_file_name, solution_file_name, solution_name, 1)

    def CollidingSpheres3D_Radiation_01(self):
        parameters_file_name = 'ProjectParametersDEM_Radiation_01.json'
        solution_name        = 'Solution_Radiation_01'
        Aux.CreateAndRunStage(KratosMultiphysics.Model(), test_folder_name, parameters_file_name, solution_file_name, solution_name, 1)

    def CollidingSpheres3D_Radiation_02(self):
        parameters_file_name = 'ProjectParametersDEM_Radiation_02.json'
        solution_name        = 'Solution_Radiation_02'
        Aux.CreateAndRunStage(KratosMultiphysics.Model(), test_folder_name, parameters_file_name, solution_file_name, solution_name, 1)

    def CollidingSpheres3D_ContactAdjust_01(self):
        parameters_file_name = 'ProjectParametersDEM_ContactAdjust_01.json'
        solution_name        = 'Solution_ContactAdjust_01'
        Aux.CreateAndRunStage(KratosMultiphysics.Model(), test_folder_name, parameters_file_name, solution_file_name, solution_name, 1)

    def CollidingSpheres3D_ContactAdjust_02(self):
        parameters_file_name = 'ProjectParametersDEM_ContactAdjust_02.json'
        solution_name        = 'Solution_ContactAdjust_02'
        Aux.CreateAndRunStage(KratosMultiphysics.Model(), test_folder_name, parameters_file_name, solution_file_name, solution_name, 1)

    def CollidingSpheres3D_ContactAdjust_03(self):
        parameters_file_name = 'ProjectParametersDEM_ContactAdjust_03.json'
        solution_name        = 'Solution_ContactAdjust_03'
        Aux.CreateAndRunStage(KratosMultiphysics.Model(), test_folder_name, parameters_file_name, solution_file_name, solution_name, 1)

if __name__ == '__main__':
    KratosUnittest.main()
