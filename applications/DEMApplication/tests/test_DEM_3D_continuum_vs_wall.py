import os
import KratosMultiphysics
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

class DEM3D_ContinuumTestVsWallSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage):

    def __init__(self, model, project_parameters):
        super().__init__(model, project_parameters)
        self.step = 0
    
    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()
        self.MoveParticle()

    def MoveParticle(self):
        for node in self.spheres_model_part.Nodes:
            V_X = V_Y = V_Z = 0.0
            if self.step == 1:
                V_X = V_Y = 0.0; V_Z = 0.15
            elif self.step == 2:
                V_X = -0.5; V_Y = V_Z = 0.0
            elif self.step == 3:
                V_X = V_Y = 0.0; V_Z = -0.15
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X, V_X)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y, V_Y)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Z, V_Z)
        self.step = self.step + 1

    @classmethod
    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM3D_continuum_vs_wall_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

class TestDEM3DContinuumVsWall(KratosUnittest.TestCase):

    def setUp(self):
        pass

    @classmethod
    def test_DEM3D_continuum_vs_wall(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM3D_continuum_vs_wall_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = KratosMultiphysics.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(DEM3D_ContinuumTestVsWallSolution, model, parameters_file_name, 1)

    def tearDown(self):
        file_to_remove = os.path.join("DEM3D_continuum_vs_wall_tests_files", "TimesPartialRelease")
        kratos_utils.DeleteFileIfExisting(GetFilePath(file_to_remove))
        os.chdir(this_working_dir_backup)

if __name__ == "__main__":
    Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
