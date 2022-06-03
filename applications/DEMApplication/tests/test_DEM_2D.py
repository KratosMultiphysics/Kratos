import os
import KratosMultiphysics as Kratos
from KratosMultiphysics import Logger
Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.DEMApplication.DEM_analysis_stage

import auxiliary_functions_for_tests

this_working_dir_backup = os.getcwd()

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class DEM2DTestSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage, KratosUnittest.TestCase):

    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM2D_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        tolerance = 1e-3
        if self.time > 0.2:
            node = self.spheres_model_part.GetNode(1)
            normal_impact_vel = node.GetSolutionStepValue(Kratos.VELOCITY_X)
            self.assertAlmostEqual(normal_impact_vel, 6.135616337653889, delta=tolerance)

            node = self.spheres_model_part.GetNode(2)
            normal_impact_vel = node.GetSolutionStepValue(Kratos.VELOCITY_X)
            self.assertAlmostEqual(normal_impact_vel, 3.532381836682557, delta=tolerance)

            node = self.spheres_model_part.GetNode(3)
            normal_impact_vel = node.GetSolutionStepValue(Kratos.VELOCITY_X)
            self.assertAlmostEqual(normal_impact_vel, 9.828777134668575, delta=tolerance)

            self.check_mark_1 = True


    def Finalize(self):
        self.assertTrue(self.check_mark_1)
        self.procedures.RemoveFoldersWithResults(str(self.main_path), str(self.problem_name), '')
        super().Finalize()

class TestDEM2D(KratosUnittest.TestCase):

    def setUp(self):
        pass

    def test_DEM2D_1(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM2D_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = Kratos.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(DEM2DTestSolution, model, parameters_file_name, auxiliary_functions_for_tests.GetHardcodedNumberOfThreads())


if __name__ == "__main__":
    Kratos.Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
