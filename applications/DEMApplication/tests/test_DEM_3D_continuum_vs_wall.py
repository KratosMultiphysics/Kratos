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

class DEM3D_ContinuumTestVsWallSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage, KratosUnittest.TestCase):

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()
        self.SetParticlePosition()

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        self.CheckForces()

    def SetParticlePosition(self):
        node = self.spheres_model_part.GetNode(1)
        node.Fix(KratosMultiphysics.VELOCITY_X)
        node.Fix(KratosMultiphysics.VELOCITY_Y)
        node.Fix(KratosMultiphysics.VELOCITY_Z)
        node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X, 0.0)
        node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y, 0.0)
        node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Z, 0.0)
        step = self.spheres_model_part.ProcessInfo.GetValue(KratosMultiphysics.TIME_STEPS)

        if step == 2:
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z, 0.2)
        elif step == 3:
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, -0.35)
        elif step == 4:
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, -0.55)
        elif step == 5:
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z, -0.3)
        elif step == 6:
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, -0.25)
        elif step == 7:
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0.0)
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z, 0.0)


        node.X = node.X0 + node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
        node.Y = node.Y0 + node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
        node.Z = node.Z0 + node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z)

    def CheckForces(self):
        node = self.spheres_model_part.GetNode(1)
        step = self.spheres_model_part.ProcessInfo.GetValue(KratosMultiphysics.TIME_STEPS)

        if step == 0:
            self.assertAlmostEqual(node.GetSolutionStepValue(DEM.CONTACT_FORCES_X), 0.0)
            self.assertAlmostEqual(node.GetSolutionStepValue(DEM.CONTACT_FORCES_Z), 0.0)
        elif step == 1:
            self.assertAlmostEqual(node.GetSolutionStepValue(DEM.CONTACT_FORCES_X), 0.0)
            self.assertAlmostEqual(node.GetSolutionStepValue(DEM.CONTACT_FORCES_Z), 0.0)
        elif step == 2:
            self.assertAlmostEqual(node.GetSolutionStepValue(DEM.CONTACT_FORCES_Z), -3065.5072065106638)
        elif step == 3:
            self.assertAlmostEqual(node.GetSolutionStepValue(DEM.CONTACT_FORCES_X), 0.0)
            self.assertAlmostEqual(node.GetSolutionStepValue(DEM.CONTACT_FORCES_Z), -3065.5072065106638)
        elif step == 4:
            self.assertAlmostEqual(node.GetSolutionStepValue(DEM.CONTACT_FORCES_X), 3065.5072065106638)
            self.assertAlmostEqual(node.GetSolutionStepValue(DEM.CONTACT_FORCES_Z), -3065.5072065106638)
        elif step == 5:
            self.assertAlmostEqual(node.GetSolutionStepValue(DEM.CONTACT_FORCES_X), 3065.5072065106638)
            self.assertAlmostEqual(node.GetSolutionStepValue(DEM.CONTACT_FORCES_Z), 0.0)
        elif step == 6:
            self.assertAlmostEqual(node.GetSolutionStepValue(DEM.CONTACT_FORCES_X), 0.0)
            self.assertAlmostEqual(node.GetSolutionStepValue(DEM.CONTACT_FORCES_Z), 0.0)
        elif step == 7:
            self.assertAlmostEqual(node.GetSolutionStepValue(DEM.CONTACT_FORCES_X), 0.0)
            self.assertAlmostEqual(node.GetSolutionStepValue(DEM.CONTACT_FORCES_Z), 0.0)

    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM3D_continuum_vs_wall_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

class TestDEM3DContinuumVsWall(KratosUnittest.TestCase):

    def setUp(self):
        pass

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
