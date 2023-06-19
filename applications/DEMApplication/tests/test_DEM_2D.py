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
        tolerance = 1e-5
        for node in self.spheres_model_part.Nodes:
            final_bouncing_vel = node.GetSolutionStepValue(Kratos.VELOCITY_X)
            if node.Id == 1:
                if self.time > 0.1999 and self.time < 0.2001:
                    self.assertAlmostEqual(final_bouncing_vel, 5.931185649037769, delta = tolerance)
            if node.Id == 2:
                if self.time > 0.1999 and self.time < 0.2001:
                    self.assertAlmostEqual(final_bouncing_vel, 7.508694851840698, delta = tolerance)
            if node.Id == 3:
                if self.time > 0.1999 and self.time < 0.2001:
                    self.assertAlmostEqual(final_bouncing_vel, 8.440679872036805, delta = tolerance)


    def Finalize(self):
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
