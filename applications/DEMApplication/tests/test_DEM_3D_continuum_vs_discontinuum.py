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

class DEM3D_ContinuumTestVsDiscontinuumSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage, KratosUnittest.TestCase):

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        self.CheckPositions()

    def CheckPositions(self):
       
        node1  = self.spheres_model_part.GetNode(1)
        node2  = self.spheres_model_part.GetNode(2)
        node3  = self.spheres_model_part.GetNode(3)
        node4  = self.spheres_model_part.GetNode(4)
        node5  = self.spheres_model_part.GetNode(5)
        node6  = self.spheres_model_part.GetNode(6)
        node7  = self.spheres_model_part.GetNode(7)
        node8  = self.spheres_model_part.GetNode(8)
        node9  = self.spheres_model_part.GetNode(9)
        node10 = self.spheres_model_part.GetNode(10)

        dipl1  = node1.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
        dipl2  = node2.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
        dipl3  = node3.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
        dipl4  = node4.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
        dipl5  = node5.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
        dipl6  = node6.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
        dipl7  = node7.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
        dipl8  = node8.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
        dipl9  = node9.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
        dipl10 = node10.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)

        self.assertAlmostEqual(dipl1, dipl3)
        self.assertAlmostEqual(dipl2, dipl4)
        self.assertAlmostEqual(dipl3, dipl5)
        self.assertAlmostEqual(dipl4, dipl6)
        self.assertAlmostEqual(dipl5, dipl7)
        self.assertAlmostEqual(dipl6, dipl8)
        self.assertAlmostEqual(dipl7, dipl9)
        self.assertAlmostEqual(dipl8, dipl10)

    @classmethod
    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM3D_continuum_vs_discontinuum_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

class TestDEM3DContinuumVsDiscontinuum(KratosUnittest.TestCase):

    def setUp(self):
        pass

    @classmethod
    def test_DEM3D_continuum_vs_discontinuum(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM3D_continuum_vs_discontinuum_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = KratosMultiphysics.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(DEM3D_ContinuumTestVsDiscontinuumSolution, model, parameters_file_name, 1)

    def tearDown(self):
        file_to_remove = os.path.join("DEM3D_continuum_vs_discontinuum_tests_files", "TimesPartialRelease")
        kratos_utils.DeleteFileIfExisting(GetFilePath(file_to_remove))
        os.chdir(this_working_dir_backup)

if __name__ == "__main__":
    Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
