import os
import KratosMultiphysics
from KratosMultiphysics import Logger
Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.DEMApplication.DEM_analysis_stage
import numpy as np

import auxiliary_functions_for_tests

this_working_dir_backup = os.getcwd()

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)


class DEM3D_SearchToleranceMain(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage, KratosUnittest.TestCase):

    def Initialize(self):
        super().Initialize()

    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_search_tolerance")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def DoCheck1(self, y_vel, tol):
        y_vel_ref = -5.8647604045771855
        self.assertAlmostEqual(y_vel, y_vel_ref, delta=tol)

    def DoCheck2(self, y_vel, tol):
        y_vel_ref = -3.3860170707636836
        self.assertAlmostEqual(y_vel, y_vel_ref, delta=tol)

    def DoCheck3(self, y_vel, tol):
        y_vel_ref = -0.5915833448135729
        self.assertAlmostEqual(y_vel, y_vel_ref, delta=tol)

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        node = self.spheres_model_part.GetNode(2)
        #reference data with freq=1 searchtolerance=0.0
        tol = 1.0e-10

        if self.spheres_model_part.ProcessInfo[KratosMultiphysics.TIME_STEPS] == 572:
            y_vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)
            self.DoCheck1(y_vel, tol)
            self.check_mark_1 = True

        if self.spheres_model_part.ProcessInfo[KratosMultiphysics.TIME_STEPS] == 3286:
            y_vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)
            self.DoCheck2(y_vel, tol)
            self.check_mark_2 = True

        if self.spheres_model_part.ProcessInfo[KratosMultiphysics.TIME_STEPS] == 6286:
            y_vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)
            self.DoCheck3(y_vel, tol)
            self.check_mark_3 = True

    def Finalize(self):
        self.assertTrue(self.check_mark_1)
        self.assertTrue(self.check_mark_2)
        self.assertTrue(self.check_mark_3)
        self.procedures.RemoveFoldersWithResults(str(self.main_path), str(self.problem_name), '')
        super().Finalize()

class DEM3D_SearchTolerance1(DEM3D_SearchToleranceMain):

    def DoCheck1(self, y_vel, tol):
        y_vel_ref = -5.866976020664773
        self.assertAlmostEqual(y_vel, y_vel_ref, delta=tol)

    def DoCheck2(self, y_vel, tol):
        y_vel_ref = -3.386926780703359
        self.assertAlmostEqual(y_vel, y_vel_ref, delta=tol)

    def DoCheck3(self, y_vel, tol):
        y_vel_ref = -0.602818980920696
        self.assertAlmostEqual(y_vel, y_vel_ref, delta=tol)

class DEM3D_SearchTolerance2(DEM3D_SearchToleranceMain):
    def DoCheck1(self, y_vel, tol):
        y_vel_ref = -5.866975875307219
        self.assertAlmostEqual(y_vel, y_vel_ref, delta=tol)

    def DoCheck2(self, y_vel, tol):
        y_vel_ref = -3.3869271295059464
        self.assertAlmostEqual(y_vel, y_vel_ref, delta=tol)

    def DoCheck3(self, y_vel, tol):
        y_vel_ref = -0.607189684524743
        self.assertAlmostEqual(y_vel, y_vel_ref, delta=tol)

class DEM3D_SearchTolerance3(DEM3D_SearchToleranceMain):

    def DoCheck1(self, y_vel, tol):
        y_vel_ref = -5.866975875307219
        self.assertAlmostEqual(y_vel, y_vel_ref, delta=tol)

    def DoCheck2(self, y_vel, tol):
        y_vel_ref = -3.38684037319764
        self.assertAlmostEqual(y_vel, y_vel_ref, delta=tol)

    def DoCheck3(self, y_vel, tol):
        y_vel_ref = -0.5971878911496257
        self.assertAlmostEqual(y_vel, y_vel_ref, delta=tol)

class TestSearchTolerance(KratosUnittest.TestCase):

    def test_SearchA(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_search_tolerance")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        with open(parameters_file_name,'r') as parameter_file:
            project_parameters = KratosMultiphysics.Parameters(parameter_file.read())
        project_parameters["SearchTolerance"].SetDouble(0.0)
        project_parameters["search_tolerance_against_walls"].SetDouble(0.0)
        project_parameters["NeighbourSearchFrequency"].SetInt(1)

        model = KratosMultiphysics.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(DEM3D_SearchToleranceMain, model, project_parameters, auxiliary_functions_for_tests.GetHardcodedNumberOfThreads())

    def test_SearchB(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_search_tolerance")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        with open(parameters_file_name,'r') as parameter_file:
            project_parameters = KratosMultiphysics.Parameters(parameter_file.read())
        project_parameters["SearchTolerance"].SetDouble(0.0)
        project_parameters["search_tolerance_against_walls"].SetDouble(0.0)
        project_parameters["NeighbourSearchFrequency"].SetInt(10)

        model = KratosMultiphysics.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(DEM3D_SearchTolerance1, model, project_parameters, auxiliary_functions_for_tests.GetHardcodedNumberOfThreads())

    def test_SearchC(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_search_tolerance")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        with open(parameters_file_name,'r') as parameter_file:
            project_parameters = KratosMultiphysics.Parameters(parameter_file.read())
        project_parameters["SearchTolerance"].SetDouble(1e-04)
        project_parameters["search_tolerance_against_walls"].SetDouble(1e-04)
        project_parameters["NeighbourSearchFrequency"].SetInt(20)

        model = KratosMultiphysics.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(DEM3D_SearchTolerance2, model, project_parameters, auxiliary_functions_for_tests.GetHardcodedNumberOfThreads())

    def test_SearchD(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_search_tolerance")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        with open(parameters_file_name,'r') as parameter_file:
            project_parameters = KratosMultiphysics.Parameters(parameter_file.read())
        project_parameters["SearchTolerance"].SetDouble(1e-03)
        project_parameters["search_tolerance_against_walls"].SetDouble(1e-03)
        project_parameters["NeighbourSearchFrequency"].SetInt(20)

        model = KratosMultiphysics.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(DEM3D_SearchTolerance3, model, project_parameters, auxiliary_functions_for_tests.GetHardcodedNumberOfThreads())


if __name__ == "__main__":
    Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
