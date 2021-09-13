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
        for node in self.spheres_model_part.Nodes:
            self.initial_normal_vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z)

    @classmethod
    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_search_tolerance")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        for node in self.spheres_model_part.Nodes:
            #reference data with freq=1 searchtolerance=0.0
            if node.Id == 2:
                tol = 1.0e-15
                if np.isclose(self.time, 0.02, rtol=0.0, atol=1e-06):
                    y_vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)
                    print(self.time, y_vel)
                    y_vel_ref = -5.86502139707038
                    self.assertAlmostEqual(y_vel, y_vel_ref, delta=tol)

                if np.isclose(self.time, 0.115, rtol=0.0, atol=1e-06):
                    y_vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)
                    print(self.time, y_vel)
                    y_vel_ref = -3.3859516373258987
                    self.assertAlmostEqual(y_vel, y_vel_ref, delta=tol)

                if np.isclose(self.time, 0.22, rtol=0.0, atol=1e-06):
                    y_vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)
                    print(self.time, y_vel)
                    y_vel_ref = -0.5929799879392164
                    self.assertAlmostEqual(y_vel, y_vel_ref, delta=tol)


    def Finalize(self):
        self.procedures.RemoveFoldersWithResults(str(self.main_path), str(self.problem_name), '')
        super().Finalize()

class DEM3D_SearchTolerance1(DEM3D_SearchToleranceMain):

    def FinalizeSolutionStep(self):
        KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage.FinalizeSolutionStep(self)
        for node in self.spheres_model_part.Nodes:
            if node.Id == 2:
                tol = 1.0e-15
                if np.isclose(self.time, 0.02, rtol=0.0, atol=1e-06):
                    y_vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)
                    print(self.time, y_vel)
                    y_vel_ref = -5.8654458179811835
                    self.assertAlmostEqual(y_vel, y_vel_ref, delta=tol)

                if np.isclose(self.time, 0.115, rtol=0.0, atol=1e-06):
                    y_vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)
                    print(self.time, y_vel)
                    y_vel_ref = -3.3861319639727263
                    self.assertAlmostEqual(y_vel, y_vel_ref, delta=tol)

                if np.isclose(self.time, 0.22, rtol=0.0, atol=1e-06):
                    y_vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)
                    print(self.time, y_vel)
                    y_vel_ref = -0.594495289987086
                    self.assertAlmostEqual(y_vel, y_vel_ref, delta=tol)

class DEM3D_SearchTolerance2(DEM3D_SearchToleranceMain):

    def FinalizeSolutionStep(self):
        KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage.FinalizeSolutionStep(self)
        for node in self.spheres_model_part.Nodes:
            if node.Id == 2:
                tol = 1.0e-15
                if np.isclose(self.time, 0.02, rtol=0.0, atol=1e-06):
                    y_vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)
                    print(self.time, y_vel)
                    y_vel_ref = -5.865445816566027
                    self.assertAlmostEqual(y_vel, y_vel_ref, delta=tol)

                if np.isclose(self.time, 0.115, rtol=0.0, atol=1e-06):
                    y_vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)
                    print(self.time, y_vel)
                    y_vel_ref = -3.386128017385994
                    self.assertAlmostEqual(y_vel, y_vel_ref, delta=tol)

                if np.isclose(self.time, 0.22, rtol=0.0, atol=1e-06):
                    y_vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)
                    print(self.time, y_vel)
                    y_vel_ref = -0.5941551772701182
                    self.assertAlmostEqual(y_vel, y_vel_ref, delta=tol)

class DEM3D_SearchTolerance3(DEM3D_SearchToleranceMain):

    def FinalizeSolutionStep(self):
        KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage.FinalizeSolutionStep(self)
        for node in self.spheres_model_part.Nodes:
            if node.Id == 2:
                tol = 1.0e-15
                if np.isclose(self.time, 0.02, rtol=0.0, atol=1e-06):
                    y_vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)
                    print(self.time, y_vel)
                    y_vel_ref = -5.86502139707038
                    self.assertAlmostEqual(y_vel, y_vel_ref, delta=tol)

                if np.isclose(self.time, 0.115, rtol=0.0, atol=1e-06):
                    y_vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)
                    print(self.time, y_vel)
                    y_vel_ref = -3.3859516373258987
                    self.assertAlmostEqual(y_vel, y_vel_ref, delta=tol)

                if np.isclose(self.time, 0.22, rtol=0.0, atol=1e-06):
                    y_vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)
                    print(self.time, y_vel)
                    y_vel_ref = -0.5929799879392164
                    self.assertAlmostEqual(y_vel, y_vel_ref, delta=tol)

class TestSearchTolerance(KratosUnittest.TestCase):

    @classmethod
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

    @classmethod
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

    @classmethod
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

    @classmethod
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
