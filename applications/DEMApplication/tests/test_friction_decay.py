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

class FrictionDecayTestSolution(DEM_analysis_stage.DEMAnalysisStage, KratosUnittest.TestCase):

    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "friction_decay_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        tolerance = 1e-4

        if self.time > 0.0499999 and self.time < 0.0500001:
            node = self.spheres_model_part.GetNode(1)
            total_force = node.GetSolutionStepValue(Kratos.TOTAL_FORCES)
            self.assertAlmostEqual(total_force[0], -4.406849, delta=tolerance)

            node = self.spheres_model_part.GetNode(2)
            total_force = node.GetSolutionStepValue(Kratos.TOTAL_FORCES)
            self.assertAlmostEqual(total_force[0], -6.480977, delta=tolerance)

            node = self.spheres_model_part.GetNode(3)
            total_force = node.GetSolutionStepValue(Kratos.TOTAL_FORCES)
            self.assertAlmostEqual(total_force[0], -4.392052, delta=tolerance)

            node = self.spheres_model_part.GetNode(4)
            total_force = node.GetSolutionStepValue(Kratos.TOTAL_FORCES)
            self.assertAlmostEqual(total_force[0], -5.724009, delta=tolerance)

            node = self.spheres_model_part.GetNode(5)
            total_force = node.GetSolutionStepValue(Kratos.TOTAL_FORCES)
            self.assertAlmostEqual(total_force[0], -4.392052, delta=tolerance)

            node = self.spheres_model_part.GetNode(6)
            total_force = node.GetSolutionStepValue(Kratos.TOTAL_FORCES)
            self.assertAlmostEqual(total_force[0], -4.406849, delta=tolerance)

            self.check_mark_1 = True


    def Finalize(self):
        self.assertTrue(self.check_mark_1)
        self.procedures.RemoveFoldersWithResults(str(self.main_path), str(self.problem_name), '')
        super().Finalize()

class TestFrictionDecay(KratosUnittest.TestCase):

    def setUp(self):
        pass

    def test_Friction_Decay(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "friction_decay_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = Kratos.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(FrictionDecayTestSolution, model, parameters_file_name, auxiliary_functions_for_tests.GetHardcodedNumberOfThreads())


if __name__ == "__main__":
    Kratos.Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
