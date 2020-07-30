import os
import KratosMultiphysics as Kratos
from KratosMultiphysics import Logger
Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
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

class DEM2DTestSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage):

    @classmethod
    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM2D_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def FinalizeTimeStep(self, time):
        tolerance = 1e-3
        for node in self.spheres_model_part.Nodes:
            normal_impact_vel = node.GetSolutionStepValue(Kratos.VELOCITY_X)
            if node.Id == 1:
                if time > 0.2:
                    expected_value = 6.076801447242313
                    self.CheckValueOfVelocity(normal_impact_vel, expected_value, tolerance)
            if node.Id == 2:
                if time > 0.2:
                    expected_value = 8.604163136887411
                    self.CheckValueOfVelocity(normal_impact_vel, expected_value, tolerance)
            if node.Id == 3:
                if time > 0.2:
                    expected_value = 10.016439272775422
                    self.CheckValueOfVelocity(normal_impact_vel, expected_value, tolerance)

    @classmethod
    def CheckValueOfVelocity(self, normal_impact_vel, expected_value, tolerance):
        if normal_impact_vel > expected_value + tolerance or normal_impact_vel < expected_value - tolerance:
            raise ValueError('Incorrect value for NORMAL_IMPACT_VELOCITY: expected value was '+ str(expected_value) + ' but received ' + str(normal_impact_vel))

    def Finalize(self):
        super(DEM2DTestSolution, self).Finalize()



class TestDEM2D(KratosUnittest.TestCase):

    def setUp(self):
        pass

    @classmethod
    def test_DEM2D_1(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM2D_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = Kratos.Model()
        CreateAndRunStageInOneOpenMPThread(DEM2DTestSolution, model, parameters_file_name)


    def tearDown(self):
        file_to_remove = os.path.join("DEM2D_tests_files", "TimesPartialRelease")
        kratos_utils.DeleteFileIfExisting(GetFilePath(file_to_remove))

        os.chdir(this_working_dir_backup)


if __name__ == "__main__":
    Kratos.Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
