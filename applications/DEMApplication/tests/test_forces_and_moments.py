import os
import KratosMultiphysics as Kratos
from KratosMultiphysics import Logger
Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.DEMApplication.DEM_analysis_stage

import KratosMultiphysics.kratos_utilities as kratos_utils
import auxiliary_functions_for_tests

this_working_dir_backup = os.getcwd()

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class ForcesAndMomentsTestSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage, KratosUnittest.TestCase):

    @classmethod
    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "forces_and_moments_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        tolerance = 1e-3
        if self.time > 0.099999 and self.time < 0.100001:
            for node in self.spheres_model_part.Nodes:
                velocity = node.GetSolutionStepValue(Kratos.VELOCITY)
                angular_velocity = node.GetSolutionStepValue(Kratos.ANGULAR_VELOCITY)
                if node.Id == 1:
                    expected_value =  0.000955
                    self.CheckValueOfVelocity(velocity, 0, expected_value, tolerance)
                    expected_value = -0.980986
                    self.CheckValueOfVelocity(velocity, 1, expected_value, tolerance)
                    expected_value = -0.000312
                    self.CheckValueOfVelocity(velocity, 2, expected_value, tolerance)
                    expected_value = -0.001194
                    self.CheckValueOfVelocity(angular_velocity, 0, expected_value, tolerance)
                    expected_value =  0.000716
                    self.CheckValueOfVelocity(angular_velocity, 1, expected_value, tolerance)
                    expected_value =  0.000718
                    self.CheckValueOfVelocity(angular_velocity, 2, expected_value, tolerance)

        if self.time > 0.299999 and self.time < 0.300001:
            for node in self.rigid_face_model_part.Nodes:
                velocity = node.GetSolutionStepValue(Kratos.VELOCITY)
                angular_velocity = node.GetSolutionStepValue(Kratos.ANGULAR_VELOCITY)
                if node.Id == 4:
                    expected_value = 154.930
                    self.CheckValueOfVelocity(velocity, 0, expected_value, tolerance)
                    expected_value = 128.674
                    self.CheckValueOfVelocity(velocity, 1, expected_value, tolerance)
                    expected_value =  21.112
                    self.CheckValueOfVelocity(velocity, 2, expected_value, tolerance)

    def CheckValueOfVelocity(self, velocity, component, expected_value, tolerance):
        self.assertAlmostEqual(velocity[component], expected_value, delta=tolerance)

    def CheckValueOfAngularVelocity(self, angular_velocity, component, expected_value, tolerance):
        self.assertAlmostEqual(angular_velocity[component], expected_value, delta=tolerance)

class TestExternalForcesAndMoments(KratosUnittest.TestCase):

    def setUp(self):
        pass

    @classmethod
    def test_ForcesAndMoments(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "forces_and_moments_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = Kratos.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(ForcesAndMomentsTestSolution, model, parameters_file_name, 1)

    def tearDown(self):
        file_to_remove = os.path.join("forces_and_moments_tests_files", "TimesPartialRelease")
        kratos_utils.DeleteFileIfExisting(GetFilePath(file_to_remove))
        file_to_remove = os.path.join("forces_and_moments_tests_files", "flux_data_new.hdf5")
        kratos_utils.DeleteFileIfExisting(GetFilePath(file_to_remove))
        os.chdir(this_working_dir_backup)


if __name__ == "__main__":
    Kratos.Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
