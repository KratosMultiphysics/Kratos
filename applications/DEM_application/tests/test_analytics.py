import os
import KratosMultiphysics as Kratos
from KratosMultiphysics import Logger
Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
from KratosMultiphysics.DEMApplication import *
import KratosMultiphysics.KratosUnittest as KratosUnittest
import main_script

import KratosMultiphysics.kratos_utilities as kratos_utils

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)


class AnalyticsTestSolution(main_script.Solution):

    def GetInputParameters(self):
        input_parameters = Kratos.Parameters(("""
        {
            "problem_name":"analytics_test_1",
            "PostNormalImpactVelocity"              : true,
            "PostTangentialImpactVelocity"          : true,
            "PostFaceNormalImpactVelocity"          : true,
            "PostFaceTangentialImpactVelocity"      : true,
            "FinalTime"                             : 0.6,
            "OutputTimeStep"                        : 1e-2
        }
        """))
        return input_parameters

    @classmethod
    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "analytics_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def FinalizeTimeStep(self, time):
        tolerance = 1e-3
        for node in self.spheres_model_part.Nodes:
            normal_impact_vel = node.GetSolutionStepValue(NORMAL_IMPACT_VELOCITY)
            face_normal_impact_vel = node.GetSolutionStepValue(FACE_NORMAL_IMPACT_VELOCITY)
            if node.Id == 1:
                if time < 0.03:
                    expected_value = 0.0
                    self.CheckValueOfNormalImpactVelocity(normal_impact_vel, expected_value, tolerance)
                    expected_value = 0.0
                    self.CheckValueOfFaceNormalImpactVelocity(face_normal_impact_vel, expected_value, tolerance)
                elif time > 0.04 and time < 0.28:
                    expected_value = 3.0
                    self.CheckValueOfNormalImpactVelocity(normal_impact_vel, expected_value, tolerance)
                    expected_value = 0.0
                    self.CheckValueOfFaceNormalImpactVelocity(face_normal_impact_vel, expected_value, tolerance)
                elif time > 0.29 and time < 0.41:
                    expected_value = 3.0
                    self.CheckValueOfNormalImpactVelocity(normal_impact_vel, expected_value, tolerance)
                    expected_value = 2.81939
                    self.CheckValueOfFaceNormalImpactVelocity(face_normal_impact_vel, expected_value, tolerance)
                elif time > 0.43:
                    expected_value = 10.268
                    self.CheckValueOfNormalImpactVelocity(normal_impact_vel, expected_value, tolerance)
                    expected_value = 7.1603
                    self.CheckValueOfFaceNormalImpactVelocity(face_normal_impact_vel, expected_value, tolerance)
            if node.Id == 2:
                if time < 0.03:
                    expected_value = 0.0
                    self.CheckValueOfNormalImpactVelocity(normal_impact_vel, expected_value, tolerance)
                elif time > 0.04 and time < 0.13:
                    expected_value = 3.0
                    self.CheckValueOfNormalImpactVelocity(normal_impact_vel, expected_value, tolerance)
                elif time > 0.15 and time < 0.42:
                    expected_value = 3.9602
                    self.CheckValueOfNormalImpactVelocity(normal_impact_vel, expected_value, tolerance)
                elif time > 0.43:
                    expected_value = 10.268
                    self.CheckValueOfNormalImpactVelocity(normal_impact_vel, expected_value, tolerance)
            if node.Id == 3:
                if time < 0.13:
                    expected_value = 0.0
                    self.CheckValueOfNormalImpactVelocity(normal_impact_vel, expected_value, tolerance)
                elif time > 0.15 and time < 0.32:
                    expected_value = 3.9602
                    self.CheckValueOfNormalImpactVelocity(normal_impact_vel, expected_value, tolerance)

    @classmethod
    def CheckValueOfNormalImpactVelocity(self, normal_impact_vel, expected_value, tolerance):
        if normal_impact_vel > expected_value + tolerance or normal_impact_vel < expected_value - tolerance:
            raise ValueError('Incorrect value for NORMAL_IMPACT_VELOCITY')

    @classmethod
    def CheckValueOfFaceNormalImpactVelocity(self, face_normal_impact_vel, expected_value, tolerance):
        if face_normal_impact_vel > expected_value + tolerance or face_normal_impact_vel < expected_value - tolerance:
            raise ValueError('Incorrect value for FACE_NORMAL_IMPACT_VELOCITY')

    def Finalize(self):
        super(AnalyticsTestSolution, self).Finalize()
        self.procedures.RemoveFoldersWithResults(self.main_path, self.problem_name)

class TestAnalytics(KratosUnittest.TestCase):

    def setUp(self):
        pass

    @classmethod
    def test_Analytics_1(self):
        AnalyticsTestSolution().Run()

    @classmethod
    def test_Analytics_2(self):
        pass

    def tearDown(self):
        file_to_remove = os.path.join("analytics_tests_files", "TimesPartialRelease")
        kratos_utils.DeleteFileIfExisting(GetFilePath(file_to_remove))


if __name__ == "__main__":
    Kratos.Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
