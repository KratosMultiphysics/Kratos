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

    @classmethod
    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "friction_decay_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        tolerance = 1e-4
        for node in self.spheres_model_part.Nodes:
            total_force = node.GetSolutionStepValue(Kratos.TOTAL_FORCES)
            if self.time > 0.0499999 and self.time < 0.0500001:
                if node.Id == 1:
                    self.assertAlmostEqual(total_force[0], -4.406849, delta=tolerance)
                if node.Id == 2:
                    self.assertAlmostEqual(total_force[0], -6.480977, delta=tolerance)
                if node.Id == 3:
                    self.assertAlmostEqual(total_force[0], -4.392052, delta=tolerance)
                if node.Id == 4:
                    self.assertAlmostEqual(total_force[0], -5.724009, delta=tolerance)
                if node.Id == 5:
                    self.assertAlmostEqual(total_force[0], -4.392052, delta=tolerance)
                if node.Id == 6:
                    self.assertAlmostEqual(total_force[0], -4.406849, delta=tolerance)


    def Finalize(self):
        self.procedures.RemoveFoldersWithResults(str(self.main_path), str(self.problem_name), '')
        super().Finalize()

class TestFrictionDecay(KratosUnittest.TestCase):

    def setUp(self):
        pass

    @classmethod
    def test_Friction_Decay(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "friction_decay_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = Kratos.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(FrictionDecayTestSolution, model, parameters_file_name, auxiliary_functions_for_tests.GetHardcodedNumberOfThreads())


if __name__ == "__main__":
    Kratos.Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
