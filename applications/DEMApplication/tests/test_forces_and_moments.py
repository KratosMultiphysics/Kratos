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

class ForcesAndMomentsTestSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage, KratosUnittest.TestCase):


    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "forces_and_moments_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def CheckValue(self, velocity, component, expected_value, tolerance):
        self.assertAlmostEqual(velocity[component], expected_value, delta=tolerance)

    def Finalize(self):
        tolerance = 1e-3
        for node in self.spheres_model_part.Nodes:
            if node.Id == 1:
                velocity = node.GetSolutionStepValue(Kratos.VELOCITY)
                expected_value = 0.000955
                self.CheckValue(velocity, 0, expected_value, tolerance)
                expected_value = 0.0
                self.CheckValue(velocity, 1, expected_value, tolerance)
                expected_value = 0.0031926545246211578
                self.CheckValue(velocity, 2, expected_value, tolerance)

            if node.Id == 2:
                angular_velocity = node.GetSolutionStepValue(Kratos.ANGULAR_VELOCITY)
                expected_value = 0.00898827541111481
                self.CheckValue(angular_velocity, 0, expected_value, tolerance)
                expected_value = 0.004494137705557405
                self.CheckValue(angular_velocity, 1, expected_value, tolerance)
                expected_value = 0.0017976550822229615
                self.CheckValue(angular_velocity, 2, expected_value, tolerance)


        for node in self.rigid_face_model_part.Nodes:
            velocity = node.GetSolutionStepValue(Kratos.VELOCITY)
            if node.Id == 4:
                expected_value = 28.810804559780266
                self.CheckValue(velocity, 0, expected_value, tolerance)
                expected_value = -7.38496
                self.CheckValue(velocity, 1, expected_value, tolerance)
                expected_value = -228.546
                self.CheckValue(velocity, 2, expected_value, tolerance)
            if node.Id == 13:
                expected_value = -4.50175
                self.CheckValue(velocity, 0, expected_value, tolerance)
                expected_value = 33.2092
                self.CheckValue(velocity, 1, expected_value, tolerance)
                expected_value = 21.0781
                self.CheckValue(velocity, 2, expected_value, tolerance)
        self.procedures.RemoveFoldersWithResults(str(self.main_path), str(self.problem_name), '')
        super().Finalize()

class TestExternalForcesAndMoments(KratosUnittest.TestCase):

    def setUp(self):
        pass


    def test_ForcesAndMoments(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "forces_and_moments_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = Kratos.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(ForcesAndMomentsTestSolution, model, parameters_file_name, auxiliary_functions_for_tests.GetHardcodedNumberOfThreads())


if __name__ == "__main__":
    Kratos.Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
