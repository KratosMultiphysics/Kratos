import os
import KratosMultiphysics
from KratosMultiphysics import Logger
Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.DEMApplication.DEM_analysis_stage

import KratosMultiphysics.kratos_utilities as kratos_utils

import auxiliary_functions_for_tests

this_working_dir_backup = os.getcwd()

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class DEM2D_ContinuumTestVsDiscontinuumSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage, KratosUnittest.TestCase):

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

        self.tolerance = 2e-4

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        self.CheckPositions()

    def CheckPositions(self):

        node1   = self.spheres_model_part.GetNode(1)
        node2   = self.spheres_model_part.GetNode(2)
        node3   = self.spheres_model_part.GetNode(3)
        node4   = self.spheres_model_part.GetNode(4)
        node5   = self.spheres_model_part.GetNode(5)
        node6   = self.spheres_model_part.GetNode(6)
        node7   = self.spheres_model_part.GetNode(7)
        node8   = self.spheres_model_part.GetNode(8)

        displ1   =  node1.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
        displ2   =  node2.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
        displ3   =  node3.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
        displ4   =  node4.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
        displ5   =  node5.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
        displ6   =  node6.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
        displ7   =  node7.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
        displ8   =  node8.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)

        self.assertAlmostEqual(displ1,   displ3, delta = self.tolerance)
        self.assertAlmostEqual(displ2,   displ4, delta = self.tolerance)
        self.assertAlmostEqual(displ5,   displ7, delta = self.tolerance)
        self.assertAlmostEqual(displ6,   displ8, delta = self.tolerance)

    @classmethod
    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM2D_continuum_vs_discontinuum_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def Finalize(self):
        self.procedures.RemoveFoldersWithResults(str(self.main_path), str(self.problem_name), '')
        super().Finalize()

class TestDEM2DContinuumVsDiscontinuum(KratosUnittest.TestCase):

    def setUp(self):
        pass

    @classmethod
    def test_DEM2D_continuum_vs_discontinuum(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM2D_continuum_vs_discontinuum_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = KratosMultiphysics.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(DEM2D_ContinuumTestVsDiscontinuumSolution, model, parameters_file_name, 1)

    def tearDown(self):
        file_to_remove = os.path.join("DEM2D_continuum_vs_discontinuum_tests_files", "TimesPartialRelease")
        kratos_utils.DeleteFileIfExisting(GetFilePath(file_to_remove))
        os.chdir(this_working_dir_backup)

if __name__ == "__main__":
    Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
