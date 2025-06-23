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

class ParallelBondModelTestSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage, KratosUnittest.TestCase):

    @classmethod
    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "dem_3d_parallel_bond_model_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        tolerance = 1e-6
        for node in self.spheres_model_part.Nodes:
            velocity = node.GetSolutionStepValue(Kratos.VELOCITY)
            angular_velocity = node.GetSolutionStepValue(Kratos.ANGULAR_VELOCITY)
            if node.Id == 1:
                if self.time > 0.0098 and self.time < 0.00981:
                    expected_value = -3.387935641635233e-05
                    self.CheckValues(velocity, 0, expected_value, tolerance)
                    expected_value = 0.0
                    self.CheckValues(velocity, 1, expected_value, tolerance)
                    expected_value = 3.20846623155654e-05
                    self.CheckValues(velocity, 2, expected_value, tolerance)
                    expected_value = 0.0
                    self.CheckValues(angular_velocity, 0, expected_value, tolerance)
                    expected_value = 7.160580099884259e-05
                    self.CheckValues(angular_velocity, 1, expected_value, tolerance)
                    expected_value = 0.0
                    self.CheckValues(angular_velocity, 2, expected_value, tolerance)

    def CheckValues(self, variable, component, expected_value, tolerance):
        self.assertAlmostEqual(variable[component], expected_value, delta=tolerance)

    def Finalize(self):
        self.procedures.RemoveFoldersWithResults(str(self.main_path), str(self.problem_name), '')
        super().Finalize()

class TestParallelBondModel(KratosUnittest.TestCase):

    def setUp(self):
        pass

    @classmethod
    def test_ParallelBondModel_1(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "dem_3d_parallel_bond_model_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = Kratos.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(ParallelBondModelTestSolution, model, parameters_file_name, 1)

if __name__ == "__main__":
    Kratos.Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
