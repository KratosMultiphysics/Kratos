import os
import KratosMultiphysics as Kratos
from KratosMultiphysics import Logger
Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.DEMApplication.DEM_analysis_stage
import numpy as np

import auxiliary_functions_for_tests

this_working_dir_backup = os.getcwd()

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class KinematicConstraintsTestSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage, KratosUnittest.TestCase):


    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "kinematic_constraints_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        tolerance = 1e-3
        for node in self.spheres_model_part.Nodes:
            velocity = node.GetSolutionStepValue(Kratos.VELOCITY)
            angular_velocity = node.GetSolutionStepValue(Kratos.ANGULAR_VELOCITY)
            if node.Id == 1:
                if self.time > 0.19 and self.time < 0.2:
                    expected_value = 0.0
                    self.CheckKinematicValues(velocity, 0, expected_value, tolerance)
                    expected_value = 0.0
                    self.CheckKinematicValues(velocity, 1, expected_value, tolerance)
                    expected_value = 0.0
                    self.CheckKinematicValues(velocity, 2, expected_value, tolerance)
                    expected_value = 0.0
                    self.CheckKinematicValues(angular_velocity, 0, expected_value, tolerance)
                    expected_value = 0.0
                    self.CheckKinematicValues(angular_velocity, 1, expected_value, tolerance)
                    expected_value = 0.0
                    self.CheckKinematicValues(angular_velocity, 2, expected_value, tolerance)
                elif np.isclose(self.time, 0.30, rtol=0.0, atol=1e-06):
                    expected_value = -0.985905
                    self.CheckKinematicValues(velocity, 1, expected_value, tolerance)

            if node.Id == 2:
                if self.time > 0.25 and self.time < 0.3:
                    expected_value = -10.0 * self.time
                    self.CheckKinematicValues(velocity, 0, expected_value, tolerance)
                if self.time > 0.3095 and self.time < 0.31:
                    expected_value = -0.09319499999999994
                    self.CheckKinematicValues(velocity, 1, expected_value, tolerance)

            if node.Id == 3:
                if self.time < 0.1 and self.time > 0.05:
                    expected_value = -5.0
                    self.CheckKinematicValues(velocity, 0, expected_value, tolerance)
                    expected_value = 0.0
                    self.CheckKinematicValues(velocity, 1, expected_value, tolerance)
                    expected_value = 0.0
                    self.CheckKinematicValues(velocity, 2, expected_value, tolerance)
                    expected_value = 0.0
                    self.CheckKinematicValues(angular_velocity, 0, expected_value, tolerance)
                    expected_value = 0.0
                    self.CheckKinematicValues(angular_velocity, 1, expected_value, tolerance)
                    expected_value = -10.0
                    self.CheckKinematicValues(angular_velocity, 2, expected_value, tolerance)

            if node.Id == 4:
                if np.isclose(self.time, 0.25, rtol=0.0, atol=1e-06):
                    expected_value = 0.21215864699077178
                    self.CheckKinematicValues(angular_velocity, 2, expected_value, tolerance)

            if node.Id == 5:
                if np.isclose(self.time, 0.30, rtol=0.0, atol=1e-06):
                    expected_value = 3.0
                    self.CheckKinematicValues(velocity, 0, expected_value, tolerance)
                    expected_value = 2.995
                    self.CheckKinematicValues(velocity, 1, expected_value, tolerance)
                    expected_value = 3.0
                    self.CheckKinematicValues(velocity, 2, expected_value, tolerance)
                    expected_value = 3.0
                    self.CheckKinematicValues(angular_velocity, 0, expected_value, tolerance)
                    expected_value = 3.0
                    self.CheckKinematicValues(angular_velocity, 1, expected_value, tolerance)
                    expected_value = 3.0
                    self.CheckKinematicValues(angular_velocity, 2, expected_value, tolerance)



    def CheckKinematicValues(self, kinematic_variable, component, expected_value, tolerance):
        self.assertAlmostEqual(kinematic_variable[component], expected_value, delta=tolerance)

    def Finalize(self):
        tolerance = 1e-3
        for node in self.rigid_face_model_part.Nodes:
            velocity = node.GetSolutionStepValue(Kratos.VELOCITY)
            angular_velocity = node.GetSolutionStepValue(Kratos.ANGULAR_VELOCITY)
            coord_x = node.X
            coord_y = node.Y
            coord_z = node.Z

            if node.Id == 10:
                    expected_value = 8.295000
                    # analytic: translation_x("1.0*t"): 8.25 + 0.045.
                    self.assertAlmostEqual(coord_x, expected_value, delta=tolerance)
                    expected_value = 0.009000
                    # analytic: initial position + rotation_x(const) + translation_z("-1.0*t**2"): - 0.25 + 0.25 + 0.009
                    self.assertAlmostEqual(coord_z, expected_value, delta=tolerance)
                    expected_value = 0.503553
                    # analytic: initial position + rotation_x(const) + translation_y(const): 0.25 + 0.103553 + 0.15
                    self.assertAlmostEqual(coord_y, expected_value, delta=tolerance)

                    expected_value = 0.30000
                    self.CheckKinematicValues(velocity, 0, expected_value, tolerance)
                    expected_value = 0.50000
                    self.CheckKinematicValues(velocity, 1, expected_value, tolerance)
                    expected_value = 1.01560
                    # analytic: v_z(t) + rotation_x(const): 0.09+2.61799*0.25*sqrt(2)
                    self.CheckKinematicValues(velocity, 2, expected_value, tolerance)
                    expected_value = 2.61799
                    self.CheckKinematicValues(angular_velocity, 0, expected_value, tolerance)

            if node.Id == 20:
                    expected_value = 2.42426
                    self.CheckKinematicValues(velocity, 0, expected_value, tolerance)
                    expected_value = -0.0523
                    self.CheckKinematicValues(velocity, 1, expected_value, tolerance)
                    expected_value = 0.65122
                    self.CheckKinematicValues(velocity, 2, expected_value, tolerance)
                    expected_value = 0.22500
                    self.CheckKinematicValues(angular_velocity, 0, expected_value, tolerance)
                    expected_value = 1.20000
                    self.CheckKinematicValues(angular_velocity, 1, expected_value, tolerance)
                    expected_value = 2.70000
                    self.CheckKinematicValues(angular_velocity, 2, expected_value, tolerance)
        self.procedures.RemoveFoldersWithResults(str(self.main_path), str(self.problem_name), '')
        super().Finalize()

class TestKinematicConstraints(KratosUnittest.TestCase):

    def setUp(self):
        pass


    def test_KinematicConstraints_1(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "kinematic_constraints_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = Kratos.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(KinematicConstraintsTestSolution, model, parameters_file_name, 1)

if __name__ == "__main__":
    Kratos.Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
