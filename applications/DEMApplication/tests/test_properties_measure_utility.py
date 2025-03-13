import os
import KratosMultiphysics as Kratos
from KratosMultiphysics import Logger
Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.DEMApplication.DEM_analysis_stage

import auxiliary_functions_for_tests

this_working_dir_backup = os.getcwd()

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class PropertiesMeasureUtilityTestSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage, KratosUnittest.TestCase):

    @classmethod
    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "properties_measure_utility_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        tolerance = 1e-8

        if self.time >= 180e-8 and self.time < 181e-8:

            #averaged stress tensor for the whole packing
            stress_tensor = self.MeasureSphereForGettingGlobalStressTensor()
            mean_stress = (stress_tensor[0][0]+stress_tensor[1][1]+stress_tensor[2][2])/3
            print("mean_stress = ", mean_stress)
            expected_value_mean_stress = 1.3623035292150967e-06
            self.assertAlmostEqual(mean_stress, expected_value_mean_stress, delta=tolerance)

    def Finalize(self):
        self.procedures.RemoveFoldersWithResults(str(self.main_path), str(self.problem_name), '')
        super().Finalize()

class TestPropertiesMeasureUtility(KratosUnittest.TestCase):

    def setUp(self):
        pass

    @classmethod
    def test_PropertiesMeasureUtility(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "properties_measure_utility_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = Kratos.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(PropertiesMeasureUtilityTestSolution, model, parameters_file_name, 1)

if __name__ == "__main__":
    Kratos.Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
