import os
import KratosMultiphysics
from KratosMultiphysics import Logger
Logger.GetDefaultOutput().SetSeverity(Logger.Severity.DETAIL)
import KratosMultiphysics.DEMApplication as DEM
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.DEMApplication.DEM_analysis_stage

import auxiliary_functions_for_tests

this_working_dir_backup = os.getcwd()

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class piecewise_linear_inletTestSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage, KratosUnittest.TestCase):

    @classmethod
    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "piecewise_linear_inlet_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        pass

    def Finalize(self):
        self.procedures.RemoveFoldersWithResults(str(self.main_path), str(self.problem_name), '')
        super().Finalize()


class TestDEM3DContact(KratosUnittest.TestCase):

    def setUp(self):
        pass

    @classmethod
    def test_piecewise_linear_inlet(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "piecewise_linear_inlet_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = KratosMultiphysics.Model()

        # Test parallel computation.
        with open(parameters_file_name,'r') as parameter_file:
            project_parameters = KratosMultiphysics.Parameters(parameter_file.read())
        piecewise_linear_inletTestSolution(model, project_parameters).Run()


if __name__ == "__main__":
    Logger.GetDefaultOutput().SetSeverity(Logger.Severity.DETAIL)
    KratosUnittest.main()
