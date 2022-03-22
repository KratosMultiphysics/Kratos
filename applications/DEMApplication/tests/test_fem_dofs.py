import os
import KratosMultiphysics as Kratos
from KratosMultiphysics import Logger
Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
import KratosMultiphysics.DEMApplication as DEM
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.DEMApplication.DEM_analysis_stage
import KratosMultiphysics.kratos_utilities as kratos_utils

import auxiliary_functions_for_tests

this_working_dir_backup = os.getcwd()

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)



class FEMDOFSTestSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage, KratosUnittest.TestCase):


    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "fem_dofs_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def Finalize(self):
        tolerance = 1e-8
        for node in self.rigid_face_model_part.Nodes:
            if node.Id == 4:
                ## right node - top fem - moves upwards - fully constrained
                velocity = node.GetSolutionStepValue(Kratos.VELOCITY_Y)
                angular_velocity = node.GetSolutionStepValue(Kratos.ANGULAR_VELOCITY_Z)
                expected_value = 1.0
                self.assertAlmostEqual(velocity, expected_value, delta=tolerance)
                expected_value = 0.0
                self.assertAlmostEqual(angular_velocity, expected_value, delta=tolerance)
            if node.Id == 8:
                ## right node - lower fem - moves downwards - partially constrained (only velocity)
                velocity = node.GetSolutionStepValue(Kratos.VELOCITY_Y)
                angular_velocity = node.GetSolutionStepValue(Kratos.ANGULAR_VELOCITY_Z)
                expected_value = -9.93245788
                self.assertAlmostEqual(velocity, expected_value, delta=tolerance)
                expected_value = -0.89324578
                self.assertAlmostEqual(angular_velocity, expected_value, delta=tolerance)

        self.procedures.RemoveFoldersWithResults(str(self.main_path), str(self.problem_name), '')
        super().Finalize()



class TestFEMDofs(KratosUnittest.TestCase):

    def setUp(self):
        pass

    def test_dofs(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "fem_dofs_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = Kratos.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(FEMDOFSTestSolution, model, parameters_file_name, 1)

if __name__ == "__main__":
    Kratos.Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
