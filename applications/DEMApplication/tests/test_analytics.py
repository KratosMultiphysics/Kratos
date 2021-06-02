import os
import KratosMultiphysics as Kratos
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



class AnalyticsTestSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage, KratosUnittest.TestCase):

    @classmethod
    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "analytics_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        tolerance = 1e-3
        for node in self.spheres_model_part.Nodes:
            normal_impact_vel = node.GetSolutionStepValue(DEM.NORMAL_IMPACT_VELOCITY)
            face_normal_impact_vel = node.GetSolutionStepValue(DEM.FACE_NORMAL_IMPACT_VELOCITY)
            if node.Id == 1:
                if self.time > 0.099:
                    expected_value = 11.07179
                    self.assertAlmostEqual(normal_impact_vel, expected_value, delta=tolerance)
                    expected_value = 6.941702
                    self.assertAlmostEqual(face_normal_impact_vel, expected_value, delta=tolerance)
            if node.Id == 2:
                if self.time > 0.099:
                    expected_value = 16.29633
                    self.assertAlmostEqual(normal_impact_vel, expected_value, delta=tolerance)
            if node.Id == 3:
                if self.time > 0.099:
                    expected_value = 16.29633
                    self.assertAlmostEqual(normal_impact_vel, expected_value, delta=tolerance)

    def Finalize(self):
        self.procedures.RemoveFoldersWithResults(str(self.main_path), str(self.problem_name), '')
        super().Finalize()

class GhostsTestSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage, KratosUnittest.TestCase):

    @classmethod
    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "analytics_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def RunAnalytics(self, time, is_time_to_print=True):
        self.MakeAnalyticsMeasurements()
        if is_time_to_print:
            self.SurfacesAnalyzerClass.CreateNewFile()
            self.SurfacesAnalyzerClass.UpdateDataBases(time)
            self.CheckTotalNumberOfCrossingParticles()

        self.SurfacesAnalyzerClass.RemoveOldFile()

    def CheckTotalNumberOfCrossingParticles(self):
        import h5py

        if self.time > 0.145:
            input_data = h5py.File(self.main_path+'/flux_data.hdf5','r')
            n_accum_h5 = input_data.get('1/n_accum')
            self.assertEqual(n_accum_h5[-1], -4)

    def Finalize(self):
        self.procedures.RemoveFoldersWithResults(str(self.main_path), str(self.problem_name), '')
        super().Finalize()

class MultiGhostsTestSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage, KratosUnittest.TestCase):

    @classmethod
    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "analytics_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def RunAnalytics(self, time, is_time_to_print=True):
        self.MakeAnalyticsMeasurements()
        if is_time_to_print:  # or IsCountStep()
            self.SurfacesAnalyzerClass.CreateNewFile()
            self.SurfacesAnalyzerClass.UpdateDataBases(time)

            if sp[Kratos.IDENTIFIER] == 'DEM-wall2':
                self.CheckTotalNumberOfCrossingParticles()

            self.SurfacesAnalyzerClass.RemoveOldFile()

    def CheckTotalNumberOfCrossingParticles(self):
        import h5py

        if self.time > 1.9:
            input_data = h5py.File(self.main_path+'/flux_data.hdf5','r')
            n_accum_h5 = input_data.get('2/n_accum')
            self.assertEqual(n_accum_h5[-1], -4)


    def Finalize(self):
        self.procedures.RemoveFoldersWithResults(str(self.main_path), str(self.problem_name), '')
        super().Finalize()

class TestAnalytics(KratosUnittest.TestCase):

    def setUp(self):
        pass

    @classmethod
    def test_Analytics_1(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "analytics_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = Kratos.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(AnalyticsTestSolution, model, parameters_file_name, 1)


    # @classmethod
    # @KratosUnittest.expectedFailure
    # def test_Analytics_2(self):
    #     path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "analytics_tests_files")
    #     parameters_file_name = os.path.join(path, "ProjectParametersDEM_single_layer_ghost.json")
    #     model = Kratos.Model()
    #     CreateAndRunStageInSelectedNumberOfOpenMPThreads(GhostsTestSolution, model, parameters_file_name, 1)


    # @classmethod
    # @KratosUnittest.expectedFailure
    # def test_Analytics_3(self):
    #     path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "analytics_tests_files")
    #     parameters_file_name = os.path.join(path, "ProjectParametersDEM_multi_layer_ghost.json")
    #      = os.path.join(path, "MaterialsDEM.json")
    #     model = Kratos.Model()
    #     CreateAndRunStageInSelectedNumberOfOpenMPThreads(MultiGhostsTestSolution, model, parameters_file_name, 1)

    def tearDown(self):
        file_to_remove = os.path.join("analytics_tests_files", "flux_data_new.hdf5")
        kratos_utils.DeleteFileIfExisting(GetFilePath(file_to_remove))
        os.chdir(this_working_dir_backup)

if __name__ == "__main__":
    Kratos.Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
