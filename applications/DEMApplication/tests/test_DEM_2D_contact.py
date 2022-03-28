import os
import KratosMultiphysics
from KratosMultiphysics import Logger
Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
import KratosMultiphysics.DEMApplication as DEM
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.DEMApplication.DEM_analysis_stage

import auxiliary_functions_for_tests

this_working_dir_backup = os.getcwd()

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class DEM2D_ContactTestSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage, KratosUnittest.TestCase):


    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM2D_contact_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def Finalize(self):
        tolerance = 1e-3
        for node in self.rigid_face_model_part.Nodes:
            dem_pressure = node.GetSolutionStepValue(DEM.DEM_PRESSURE)
            contact_force = node.GetSolutionStepValue(DEM.CONTACT_FORCES_Y)
            if node.Id == 13:
                expected_value = 45072.567
                self.assertAlmostEqual(dem_pressure, expected_value, delta=tolerance)
                expected_value = -23114.137
                self.assertAlmostEqual(contact_force, expected_value, delta=tolerance)
            if node.Id == 22:
                expected_value = 26774.109
                self.assertAlmostEqual(dem_pressure, expected_value, delta=tolerance)
                expected_value = -13730.312
                self.assertAlmostEqual(contact_force, expected_value, delta=tolerance)
        self.procedures.RemoveFoldersWithResults(str(self.main_path), str(self.problem_name), '')
        super().Finalize()


class TestDEM2DContact(KratosUnittest.TestCase):

    def setUp(self):
        pass


    def test_DEM2D_contact(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM2D_contact_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = KratosMultiphysics.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(DEM2D_ContactTestSolution, model, parameters_file_name, auxiliary_functions_for_tests.GetHardcodedNumberOfThreads())

if __name__ == "__main__":
    Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
