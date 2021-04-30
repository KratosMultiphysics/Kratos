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

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        self.SetParticlePosition()

    def SetParticlePosition(self):
        for node in self.spheres_model_part.Nodes:
            if self.step == 0:
                node.X = 0.75
                node.Y = 0.7
                node.Z = 0.6
            elif self.step == 1:
                node.X = 0.75
                node.Y = 0.7
                node.Z = 0.8
            elif self.step == 2:
                node.X = 0.4
                node.Y = 0.7
                node.Z = 0.8
            elif self.step == 3:
                node.X = 0.3
                node.Y = 0.7
                node.Z = 0.8
            elif self.step == 4:
                node.X = 0.3
                node.Y = 0.7
                node.Z = 0.6
            elif self.step == 5:
                node.X = 0.3
                node.Y = 0.7
                node.Z = 0.3
            elif self.step == 6:
                node.X = 0.5
                node.Y = 0.7
                node.Z = 0.3
            elif self.step == 7:
                node.X = 0.75
                node.Y = 0.7
                node.Z = 0.6

            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X, 0.0)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y, 0.0)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Z, 0.0)
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
