import os
import KratosMultiphysics
from KratosMultiphysics import Logger
Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
import KratosMultiphysics.DEMApplication as DEM
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.DEMApplication.DEM_analysis_stage import DEMAnalysisStage

import KratosMultiphysics.kratos_utilities as kratos_utils

import auxiliary_functions_for_tests

this_working_dir_backup = os.getcwd()

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class DEM2D_ControlModuleTestSolution(DEMAnalysisStage, KratosUnittest.TestCase):

    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM2D_control_module_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def Initialize(self):
        super().Initialize()

        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM2D_control_module_tests_files")
        cm_project_parameters_file_name = os.path.join(path, "cm_parameters.json")

        with open(cm_project_parameters_file_name,'r') as parameters_file:
            self.cm_project_parameters = KratosMultiphysics.Parameters(parameters_file.read())

        #NOTE: We will transform CM utility into a process eventually
        from KratosMultiphysics.DEMApplication.multiaxial_control_module_generalized_2d_utility import MultiaxialControlModuleGeneralized2DUtility
        self.multiaxial_control_module = MultiaxialControlModuleGeneralized2DUtility(self.model, self.cm_project_parameters)
        self.multiaxial_control_module.ExecuteInitialize()

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

        self.multiaxial_control_module.ExecuteInitializeSolutionStep()

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

        self.multiaxial_control_module.ExecuteFinalizeSolutionStep()

    def PrintResultsForGid(self, time):
        super().PrintResultsForGid(time)

        self.multiaxial_control_module.PrintResults()

    def Finalize(self):
        tolerance = 1.001
        for node in self.rigid_face_model_part.Nodes:
            if node.Id == 5:
                node_force_x = node.GetSolutionStepValue(DEM.CONTACT_FORCES_X)
                expected_value = 316.79
                self.assertAlmostEqual(node_force_x, expected_value, delta=tolerance)
            elif node.Id == 6:
                node_force_y = node.GetSolutionStepValue(DEM.CONTACT_FORCES_Y)
                expected_value = 150.1
                self.assertAlmostEqual(node_force_y, expected_value, delta=tolerance)

        super().Finalize()

class TestDEM2DControlModule(KratosUnittest.TestCase):

    def setUp(self):
        pass

    def test_DEM2D_control_module(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM2D_control_module_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = KratosMultiphysics.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(DEM2D_ControlModuleTestSolution, model, parameters_file_name, 1)

    def tearDown(self):
        file_to_remove = os.path.join("DEM2D_control_module_tests_files", "TimesPartialRelease")
        kratos_utils.DeleteFileIfExisting(GetFilePath(file_to_remove))
        os.chdir(this_working_dir_backup)


if __name__ == "__main__":
    Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
