import os
import KratosMultiphysics as Kratos
import KratosMultiphysics.DEMApplication as DEM
from KratosMultiphysics import Logger
Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.DEMApplication.DEM_analysis_stage

import KratosMultiphysics.kratos_utilities as kratos_utils
import auxiliary_functions_for_tests

this_working_dir_backup = os.getcwd()

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class HistoryDependentCLsTestSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage, KratosUnittest.TestCase):

    @classmethod
    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "history_dependent_CLs_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        tolerance = 1e-3
        if self.time > 0.09999 and self.time < 0.10001:
            for node in self.spheres_model_part.Nodes:
                force = node.GetSolutionStepValue(DEM.CONTACT_FORCES)

                if node.Id == 2:
                    expected_value = -5.03801
                    self.CheckValueOfForce(force, 0, expected_value, tolerance)
                    expected_value = 0.0
                    self.CheckValueOfForce(force, 1, expected_value, tolerance)
                    expected_value = 44.26305
                    self.CheckValueOfForce(force, 2, expected_value, tolerance)

                if node.Id == 3:
                    expected_value = -7.37913
                    self.CheckValueOfForce(force, 0, expected_value, tolerance)
                    expected_value = 0.0
                    self.CheckValueOfForce(force, 1, expected_value, tolerance)
                    expected_value = 55.2373
                    self.CheckValueOfForce(force, 2, expected_value, tolerance)

                if node.Id == 5:
                    expected_value = -1219.77524
                    self.CheckValueOfForce(force, 0, expected_value, tolerance)
                    expected_value = 0.0
                    self.CheckValueOfForce(force, 1, expected_value, tolerance)
                    expected_value = 1854.85877
                    self.CheckValueOfForce(force, 2, expected_value, tolerance)

                if node.Id == 6:
                    expected_value = -2646.11603
                    self.CheckValueOfForce(force, 0, expected_value, tolerance)
                    expected_value = 0.0
                    self.CheckValueOfForce(force, 1, expected_value, tolerance)
                    expected_value = 3850.74485
                    self.CheckValueOfForce(force, 2, expected_value, tolerance)

    def CheckValueOfForce(self, force, component, expected_value, tolerance):
        self.assertAlmostEqual(force[component], expected_value, delta=tolerance)

class TestHistoryDependentCLs(KratosUnittest.TestCase):

    def setUp(self):
        pass

    @classmethod
    def test_HistoryDependentCLs(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "history_dependent_CLs_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = Kratos.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(HistoryDependentCLsTestSolution, model, parameters_file_name, 1)

    def tearDown(self):
        file_to_remove = os.path.join("history_dependent_CLs_tests_files", "TimesPartialRelease")
        kratos_utils.DeleteFileIfExisting(GetFilePath(file_to_remove))
        file_to_remove = os.path.join("history_dependent_CLs_tests_files", "flux_data_new.hdf5")
        kratos_utils.DeleteFileIfExisting(GetFilePath(file_to_remove))
        os.chdir(this_working_dir_backup)


if __name__ == "__main__":
    Kratos.Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
