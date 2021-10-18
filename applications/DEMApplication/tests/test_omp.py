import os
import KratosMultiphysics
from KratosMultiphysics import Logger
Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.DEMApplication.DEM_analysis_stage

import auxiliary_functions_for_tests

this_working_dir_backup = os.getcwd()

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class OMPTestSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage, KratosUnittest.TestCase):


    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_omp")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def CheckValues(self, y_disp, y_vel):
        tol = 1.0e-13
        # DEM reference values
        y_velr  = 8.1296828952707
        y_dispr = 0.11937225925581423

        self.assertAlmostEqual(y_disp, y_dispr, delta=tol)
        self.assertAlmostEqual(y_vel, y_velr, delta=tol)

    def Finalize(self):
        for node in self.spheres_model_part.Nodes:
            if node.Id == 5:
                y_disp = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
                y_vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)
                print(y_vel)

        # self.CheckValues(y_disp, y_vel)
        Logger.PrintWarning("end thread")
        Logger.Flush()
        self.procedures.RemoveFoldersWithResults(str(self.main_path), str(self.problem_name), '')
        super().Finalize()


class TestOMP(KratosUnittest.TestCase):


    def setUp(self):
        pass

    def test_omp_1(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_omp")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = KratosMultiphysics.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(OMPTestSolution, model, parameters_file_name, 1)


    def test_omp_2(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_omp")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = KratosMultiphysics.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(OMPTestSolution, model, parameters_file_name, 2)


    def test_omp_3(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_omp")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = KratosMultiphysics.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(OMPTestSolution, model, parameters_file_name, 3)


    def test_omp_4(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_omp")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = KratosMultiphysics.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(OMPTestSolution, model, parameters_file_name, 4)


if __name__ == "__main__":
    Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
