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

class DEM3D_ContinuumTestVsDiscontinuumSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage, KratosUnittest.TestCase):

    def Initialize(self):
        super().Initialize()
        self._GetSolver().cplusplus_strategy.BreakAllBonds()
    
    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

        self.tolerance = 5e-5

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
        node9   = self.spheres_model_part.GetNode(9)
        node10  = self.spheres_model_part.GetNode(10)
        node11  = self.spheres_model_part.GetNode(11)
        node12  = self.spheres_model_part.GetNode(12)
        node13  = self.spheres_model_part.GetNode(13)
        node14  = self.spheres_model_part.GetNode(14)
        node15  = self.spheres_model_part.GetNode(15)
        node16  = self.spheres_model_part.GetNode(16)

        displ1   =  node1.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
        displ2   =  node2.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
        displ3   =  node3.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
        displ4   =  node4.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
        displ5   =  node5.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
        displ6   =  node6.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
        displ7   =  node7.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
        displ8   =  node8.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
        displ9   =  node9.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
        displ10  = node10.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
        displ11  = node11.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
        displ12  = node12.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
        displ13  = node13.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
        displ14  = node14.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
        displ15  = node15.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
        displ16  = node16.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)

        self.assertAlmostEqual(displ1,   displ3, delta = self.tolerance)
        self.assertAlmostEqual(displ2,   displ4, delta = self.tolerance)
        self.assertAlmostEqual(displ5,   displ7, delta = self.tolerance)
        self.assertAlmostEqual(displ6,   displ8, delta = self.tolerance)
        self.assertAlmostEqual(displ9,  displ11, delta = self.tolerance)
        self.assertAlmostEqual(displ10, displ12, delta = self.tolerance)
        self.assertAlmostEqual(displ13, displ15, delta = self.tolerance)
        self.assertAlmostEqual(displ14, displ16, delta = self.tolerance)

    @classmethod
    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM3D_continuum_vs_discontinuum_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def Finalize(self):
        self.procedures.RemoveFoldersWithResults(str(self.main_path), str(self.problem_name), '')
        super().Finalize()


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
