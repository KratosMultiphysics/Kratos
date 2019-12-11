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

class KinematicConstraintsTestSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage):

    @classmethod
    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "kinematic_constraints_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def FinalizeTimeStep(self, time):
        tolerance = 1e-3
        for node in self.spheres_model_part.Nodes:
            velocity = node.GetSolutionStepValue(Kratos.VELOCITY)
            angular_velocity = node.GetSolutionStepValue(Kratos.ANGULAR_VELOCITY)
            if node.Id == 1:
                if time > 0.18 and time < 0.2:
                    expected_value = 0.0
                    self.CheckValueOfVelocity(velocity, 0, expected_value, tolerance)
                    expected_value = 0.0
                    self.CheckValueOfVelocity(velocity, 1, expected_value, tolerance)
                    expected_value = 0.0
                    self.CheckValueOfVelocity(velocity, 2, expected_value, tolerance)
                    expected_value = 0.0
                    self.CheckValueOfAngularVelocity(angular_velocity, 0, expected_value, tolerance)
                    expected_value = 0.0
                    self.CheckValueOfAngularVelocity(angular_velocity, 1, expected_value, tolerance)
                    expected_value = 0.0
                    self.CheckValueOfAngularVelocity(angular_velocity, 2, expected_value, tolerance)
                elif time > 0.31999 and time < 0.32:
                    expected_value = -1.179
                    self.CheckValueOfVelocity(velocity, 1, expected_value, tolerance)

            if node.Id == 2:
                if time > 0.25 and time < 0.3:
                    expected_value = -10.0 * time
                    self.CheckValueOfVelocity(velocity, 0, expected_value, tolerance)
                if time > 0.59999 and time < 0.6:
                    expected_value = -1.962
                    self.CheckValueOfVelocity(velocity, 1, expected_value, tolerance)

            if node.Id == 3:
                if time < 0.1:
                    expected_value = -5.0
                    self.CheckValueOfVelocity(velocity, 0, expected_value, tolerance)
                    expected_value = 0.0
                    self.CheckValueOfVelocity(velocity, 1, expected_value, tolerance)
                    expected_value = 0.0
                    self.CheckValueOfVelocity(velocity, 2, expected_value, tolerance)
                    expected_value = 0.0
                    self.CheckValueOfAngularVelocity(angular_velocity, 0, expected_value, tolerance)
                    expected_value = 0.0
                    self.CheckValueOfAngularVelocity(angular_velocity, 1, expected_value, tolerance)
                    expected_value = -10.0
                    self.CheckValueOfAngularVelocity(angular_velocity, 2, expected_value, tolerance)
            if node.Id == 4:
                if time > 0.22 and time < 0.25:
                    expected_value = 0.2192
                    self.CheckValueOfAngularVelocity(angular_velocity, 2, expected_value, tolerance)

    @classmethod
    def CheckValueOfVelocity(self, velocity, component, expected_value, tolerance):
        if velocity[component] > expected_value + tolerance or velocity[component] < expected_value - tolerance:
            raise ValueError('Incorrect value for VELOCITY ' + str(component) + ': expected value was '+ str(expected_value) + ' but received ' + str(velocity))

    @classmethod
    def CheckValueOfAngularVelocity(self, angular_velocity, component, expected_value, tolerance):
        if angular_velocity[component] > expected_value + tolerance or angular_velocity[component] < expected_value - tolerance:
            raise ValueError('Incorrect value for ANGULAR_VELOCITY ' + str(component) + ': expected value was '+ str(expected_value) + ' but received ' + str(angular_velocity))

    def Finalize(self):
        super(KinematicConstraintsTestSolution, self).Finalize()


class TestKinematicConstraints(KratosUnittest.TestCase):

    def setUp(self):
        pass

    @classmethod
    def test_KinematicConstraints_1(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "kinematic_constraints_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = Kratos.Model()
        CreateAndRunStageInOneOpenMPThread(KinematicConstraintsTestSolution, model, parameters_file_name)


    def tearDown(self):
        file_to_remove = os.path.join("kinematic_constraints_tests_files", "TimesPartialRelease")
        kratos_utils.DeleteFileIfExisting(GetFilePath(file_to_remove))

        os.chdir(this_working_dir_backup)


if __name__ == "__main__":
    Kratos.Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
