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

    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "history_dependent_CLs_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        tolerance = 1e-3
        if self.time > 0.09999 and self.time < 0.10001:
            node = self.spheres_model_part.GetNode(2)
            force = node.GetSolutionStepValue(DEM.CONTACT_FORCES)
            expected_value = -5.0380
            self.CheckValueOfForce(force, 0, expected_value, tolerance)
            expected_value = 0.0
            self.CheckValueOfForce(force, 1, expected_value, tolerance)
            expected_value = 44.2631
            self.CheckValueOfForce(force, 2, expected_value, tolerance)

            node = self.spheres_model_part.GetNode(3)
            force = node.GetSolutionStepValue(DEM.CONTACT_FORCES)
            expected_value = -7.3791
            self.CheckValueOfForce(force, 0, expected_value, tolerance)
            expected_value = 0.0
            self.CheckValueOfForce(force, 1, expected_value, tolerance)
            expected_value = 55.2373
            self.CheckValueOfForce(force, 2, expected_value, tolerance)

            node = self.spheres_model_part.GetNode(5)
            force = node.GetSolutionStepValue(DEM.CONTACT_FORCES)
            expected_value = -1152.7767
            self.CheckValueOfForce(force, 0, expected_value, tolerance)
            expected_value = 0.0
            self.CheckValueOfForce(force, 1, expected_value, tolerance)
            expected_value = 1752.9904
            self.CheckValueOfForce(force, 2, expected_value, tolerance)

            node = self.spheres_model_part.GetNode(6)
            force = node.GetSolutionStepValue(DEM.CONTACT_FORCES)
            expected_value = -2501.1437
            self.CheckValueOfForce(force, 0, expected_value, tolerance)
            expected_value = 0.0
            self.CheckValueOfForce(force, 1, expected_value, tolerance)
            expected_value = 3639.8021
            self.CheckValueOfForce(force, 2, expected_value, tolerance)

            self.check_mark_1 = True

    def CheckValueOfForce(self, force, component, expected_value, tolerance):
        self.assertAlmostEqual(force[component], expected_value, delta=tolerance)

    def Finalize(self):
        self.assertTrue(self.check_mark_1)
        super().Finalize()

class TestHistoryDependentCLs(KratosUnittest.TestCase):

    def setUp(self):
        pass

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
