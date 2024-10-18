import os
import KratosMultiphysics as Kratos
from KratosMultiphysics import Logger
Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.DEMApplication.DEM_analysis_stage as DEM_analysis_stage

import auxiliary_functions_for_tests

this_working_dir_backup = os.getcwd()

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class GluedParticlesTestSolution(DEM_analysis_stage.DEMAnalysisStage, KratosUnittest.TestCase):

    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "glued_particles_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        tolerance = 1e-4
        node = self.spheres_model_part.GetNode(1)
        angular_velocity = node.GetSolutionStepValue(Kratos.ANGULAR_VELOCITY)

        if self.time > 0.01:
            self.assertAlmostEqual(angular_velocity[0], 2.0, delta=tolerance)
            self.check_mark_1 = True

        if self.time > 0.0499999 and self.time < 0.05000001:
            self.assertAlmostEqual(node.X, -1.0, delta=tolerance)
            self.assertAlmostEqual(node.Y, 0.960067, delta=tolerance)
            self.assertAlmostEqual(node.Z, 0.398002, delta=tolerance)
            self.check_mark_2 = True

        if self.time > 0.0999999 and self.time < 0.10000001:
            self.assertAlmostEqual(node.X, -1.0, tolerance)
            self.assertAlmostEqual(node.Y, 0.920532, delta=tolerance)
            self.assertAlmostEqual(node.Z, 0.392027, delta=tolerance)
            self.check_mark_3 = True


    def Finalize(self):
        self.assertTrue(self.check_mark_1)
        self.assertTrue(self.check_mark_2)
        self.assertTrue(self.check_mark_3)
        self.procedures.RemoveFoldersWithResults(str(self.main_path), str(self.problem_name), '')
        super().Finalize()

class TestGluedParticles(KratosUnittest.TestCase):

    def setUp(self):
        pass

    def test_Glued_Particles_1(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "glued_particles_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = Kratos.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(GluedParticlesTestSolution, model, parameters_file_name, auxiliary_functions_for_tests.GetHardcodedNumberOfThreads())


if __name__ == "__main__":
    Kratos.Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
