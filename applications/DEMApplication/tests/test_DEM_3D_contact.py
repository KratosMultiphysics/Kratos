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

class DEM3D_ContactTestSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage, KratosUnittest.TestCase):


    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM3D_contact_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def Finalize(self):
        tolerance = 1e-3
        for node in self.rigid_face_model_part.Nodes:
            dem_pressure = node.GetSolutionStepValue(DEM.DEM_PRESSURE)
            contact_force = node.GetSolutionStepValue(DEM.CONTACT_FORCES_Z)
            if node.Id == 9:
                self.assertAlmostEqual(dem_pressure, 2850.829, delta=tolerance)
                self.assertAlmostEqual(contact_force, -11403.317, delta=tolerance)
            if node.Id == 13:
                self.assertAlmostEqual(dem_pressure, 1273.402, delta=tolerance)
                self.assertAlmostEqual(contact_force, -5093.607, delta=tolerance)
        self.procedures.RemoveFoldersWithResults(str(self.main_path), str(self.problem_name), '')
        super().Finalize()


class TestDEM3DContact(KratosUnittest.TestCase):

    def setUp(self):
        pass


    def test_DEM3D_contact(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM3D_contact_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = KratosMultiphysics.Model()

        # Test parallel computation.
        with open(parameters_file_name,'r') as parameter_file:
            project_parameters = KratosMultiphysics.Parameters(parameter_file.read())
        DEM3D_ContactTestSolution(model, project_parameters).Run()


if __name__ == "__main__":
    Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
