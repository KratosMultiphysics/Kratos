import os
import KratosMultiphysics as Kratos
from KratosMultiphysics import Logger
Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
import KratosMultiphysics.DEMApplication as DEM
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.DEMApplication.DEM_analysis_stage

import KratosMultiphysics.kratos_utilities as kratos_utils

this_working_dir_backup = os.getcwd()

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

def CreateAndRunStageInOneOpenMPThread(my_obj, model, parameters_file_name):
    omp_utils = Kratos.OpenMPUtils()
    if "OMP_NUM_THREADS" in os.environ:
        initial_number_of_threads = os.environ['OMP_NUM_THREADS']
        omp_utils.SetNumThreads(1)

    with open(parameters_file_name,'r') as parameter_file:
        project_parameters = Kratos.Parameters(parameter_file.read())

    my_obj(model, project_parameters).Run()

    if "OMP_NUM_THREADS" in os.environ:
        omp_utils.SetNumThreads(int(initial_number_of_threads))

class AnalyticsTestSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage):

    @classmethod
    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "analytics_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def FinalizeTimeStep(self, time):
        tolerance = 1e-3
        for node in self.spheres_model_part.Nodes:
            normal_impact_vel = node.GetSolutionStepValue(DEM.NORMAL_IMPACT_VELOCITY)
            face_normal_impact_vel = node.GetSolutionStepValue(DEM.FACE_NORMAL_IMPACT_VELOCITY)
            if node.Id == 1:
                if time > 0.099:
                    expected_value = 11.07179
                    self.CheckValueOfNormalImpactVelocity(normal_impact_vel, expected_value, tolerance)
                    expected_value = 6.941702
                    self.CheckValueOfFaceNormalImpactVelocity(face_normal_impact_vel, expected_value, tolerance)
            if node.Id == 2:
                if time > 0.099:
                    expected_value = 16.29633
                    self.CheckValueOfNormalImpactVelocity(normal_impact_vel, expected_value, tolerance)
            if node.Id == 3:
                if time > 0.099:
                    expected_value = 16.29633
                    self.CheckValueOfNormalImpactVelocity(normal_impact_vel, expected_value, tolerance)

    @classmethod
    def CheckValueOfNormalImpactVelocity(self, normal_impact_vel, expected_value, tolerance):
        if normal_impact_vel > expected_value + tolerance or normal_impact_vel < expected_value - tolerance:
            raise ValueError('Incorrect value for NORMAL_IMPACT_VELOCITY: expected value was '+ str(expected_value) + ' but received ' + str(normal_impact_vel))

    @classmethod
    def CheckValueOfFaceNormalImpactVelocity(self, face_normal_impact_vel, expected_value, tolerance):
        if face_normal_impact_vel > expected_value + tolerance or face_normal_impact_vel < expected_value - tolerance:
            raise ValueError('Incorrect value for FACE_NORMAL_IMPACT_VELOCITY: expected value was '+ str(expected_value) + ' but received ' + str(face_normal_impact_vel))

    def Finalize(self):
        super(AnalyticsTestSolution, self).Finalize()


class GhostsTestSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage):

    @classmethod
    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "analytics_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def RunAnalytics(self, time, is_time_to_print=True):
            self.MakeAnalyticsMeasurements()
            if is_time_to_print:
                self.FaceAnalyzerClass.CreateNewFile()
                for sp in (sp for sp in self.rigid_face_model_part.SubModelParts if sp[DEM.IS_GHOST]):
                    self.face_watcher_analysers[sp.Name].UpdateDataFiles(time)
                    self.CheckTotalNumberOfCrossingParticles()

            self.FaceAnalyzerClass.RemoveOldFile()

    def CheckTotalNumberOfCrossingParticles(self):
        import h5py

        if self.time > 0.145:
            input_data = h5py.File(self.main_path+'/flux_data.hdf5','r')
            n_accum_h5 = input_data.get('1/n_accum')

            if n_accum_h5[-1] != -4:
                print(n_accum_h5[-1])
                raise ValueError('The total value of crossing particles was not the expected!')

    def Finalize(self):
        super(GhostsTestSolution, self).Finalize()



class MultiGhostsTestSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage):

    @classmethod
    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "analytics_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def RunAnalytics(self, time, is_time_to_print=True):
            self.MakeAnalyticsMeasurements()
            if is_time_to_print:  # or IsCountStep()
                self.FaceAnalyzerClass.CreateNewFile()
                for sp in (sp for sp in self.rigid_face_model_part.SubModelParts if sp[DEM.IS_GHOST]):
                    self.face_watcher_analysers[sp.Name].UpdateDataFiles(time)

                    if sp[Kratos.IDENTIFIER] == 'DEM-wall2':
                        self.CheckTotalNumberOfCrossingParticles()

                self.FaceAnalyzerClass.RemoveOldFile()

    def CheckTotalNumberOfCrossingParticles(self):
        import h5py

        if self.time > 1.9:
            input_data = h5py.File(self.main_path+'/flux_data.hdf5','r')
            n_accum_h5 = input_data.get('2/n_accum')

            if n_accum_h5[-1] != -4:
                print(n_accum_h5[-1])
                raise ValueError('The total value of crossing particles was not the expected!')

    def Finalize(self):
        super(MultiGhostsTestSolution, self).Finalize()


class TestAnalytics(KratosUnittest.TestCase):

    def setUp(self):
        pass

    @classmethod
    def test_Analytics_1(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "analytics_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = Kratos.Model()
        CreateAndRunStageInOneOpenMPThread(AnalyticsTestSolution, model, parameters_file_name)


    # @classmethod
    # @KratosUnittest.expectedFailure
    # def test_Analytics_2(self):
    #     path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "analytics_tests_files")
    #     parameters_file_name = os.path.join(path, "ProjectParametersDEM_single_layer_ghost.json")
    #     model = Kratos.Model()
    #     CreateAndRunStageInOneOpenMPThread(GhostsTestSolution, model, parameters_file_name)


    # @classmethod
    # @KratosUnittest.expectedFailure
    # def test_Analytics_3(self):
    #     path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "analytics_tests_files")
    #     parameters_file_name = os.path.join(path, "ProjectParametersDEM_multi_layer_ghost.json")
    #     model = Kratos.Model()
    #     CreateAndRunStageInOneOpenMPThread(MultiGhostsTestSolution, model, parameters_file_name)


    def tearDown(self):
        file_to_remove = os.path.join("analytics_tests_files", "TimesPartialRelease")
        kratos_utils.DeleteFileIfExisting(GetFilePath(file_to_remove))

        os.chdir(this_working_dir_backup)


if __name__ == "__main__":
    Kratos.Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
