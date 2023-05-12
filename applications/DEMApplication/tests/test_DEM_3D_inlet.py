import os
import KratosMultiphysics
from KratosMultiphysics import Logger
Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
import KratosMultiphysics.DEMApplication as DEM
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.DEMApplication.DEM_analysis_stage

import auxiliary_functions_for_tests

this_working_dir_backup = os.getcwd()

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class DEM3D_InletTestSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage, KratosUnittest.TestCase):

    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM3D_inlet_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        tolerance = 1e-8

        if self.time >= 1.15:
            node = self.spheres_model_part.GetNode(10)
            node_vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)
            node_force = node.GetSolutionStepValue(KratosMultiphysics.TOTAL_FORCES_Y)
            self.assertAlmostEqual(node_vel, 0.40632795240200154, delta=tolerance)
            self.assertAlmostEqual(node_force, -68713.6675439023, delta=tolerance)

            node = self.spheres_model_part.GetNode(11)
            node_disp = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
            self.assertAlmostEqual(node_disp, 0.046134545389347956, delta=tolerance)

            self.check_mark_1 = True


    def Finalize(self):
        self.assertTrue(self.check_mark_1)
        self.procedures.RemoveFoldersWithResults(str(self.main_path), str(self.problem_name), '')
        super().Finalize()


class TestDEM3DInlet(KratosUnittest.TestCase):

    def setUp(self):
        pass

    def test_DEM3D_inlet(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM3D_inlet_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = KratosMultiphysics.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(DEM3D_InletTestSolution, model, parameters_file_name, auxiliary_functions_for_tests.GetHardcodedNumberOfThreads())


if __name__ == "__main__":
    Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
