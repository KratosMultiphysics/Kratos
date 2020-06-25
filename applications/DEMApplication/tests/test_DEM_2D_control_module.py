import os
import KratosMultiphysics
from KratosMultiphysics import Logger
Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
import KratosMultiphysics.DEMApplication as DEM
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.DEMApplication.DEM_analysis_stage import DEMAnalysisStage

import KratosMultiphysics.kratos_utilities as kratos_utils

this_working_dir_backup = os.getcwd()

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

def CreateAndRunStageInOneOpenMPThread(my_obj, model, parameters_file_name):

    # os.environ['OMP_NUM_THREADS']='1'
    # os.system("echo Unittest will be running on $OMP_NUM_THREADS cpu")

    omp_utils = KratosMultiphysics.OpenMPUtils()
    if "OMP_NUM_THREADS" in os.environ:
        initial_number_of_threads = os.environ['OMP_NUM_THREADS']
        omp_utils.SetNumThreads(1)

    with open(parameters_file_name,'r') as parameter_file:
        project_parameters = KratosMultiphysics.Parameters(parameter_file.read())

    my_obj(model, project_parameters).Run()

    if "OMP_NUM_THREADS" in os.environ:
        omp_utils.SetNumThreads(int(initial_number_of_threads))



class DEM2D_ControlModuleTestSolution(DEMAnalysisStage):

    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM2D_control_module_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def Initialize(self):
        super(DEM2D_ControlModuleTestSolution, self).Initialize()

        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM2D_control_module_tests_files")
        cm_project_parameters_file_name = os.path.join(path, "cm_parameters.json")

        with open(cm_project_parameters_file_name,'r') as parameters_file:
            self.cm_project_parameters = KratosMultiphysics.Parameters(parameters_file.read())

        #NOTE: We will transform CM utility into a process eventually
        from KratosMultiphysics.DEMApplication.multiaxial_control_module_generalized_2d_utility import MultiaxialControlModuleGeneralized2DUtility
        self.multiaxial_control_module = MultiaxialControlModuleGeneralized2DUtility(self.model, self.cm_project_parameters)
        self.multiaxial_control_module.ExecuteInitialize()

    def InitializeSolutionStep(self):
        super(DEM2D_ControlModuleTestSolution, self).InitializeSolutionStep()

        self.multiaxial_control_module.ExecuteInitializeSolutionStep()

    def FinalizeSolutionStep(self):
        super(DEM2D_ControlModuleTestSolution, self).FinalizeSolutionStep()

        self.multiaxial_control_module.ExecuteFinalizeSolutionStep()

    def PrintResultsForGid(self, time):
        super(DEM2D_ControlModuleTestSolution, self).PrintResultsForGid(time)

        self.multiaxial_control_module.PrintResults()

    def Finalize(self):
        tolerance = 1.001
        for node in self.rigid_face_model_part.Nodes:
            if node.Id == 5:
                node_force_x = node.GetSolutionStepValue(DEM.CONTACT_FORCES_X)
                expected_value = 316.79
                self.CheckForceX(node_force_x, expected_value, tolerance)
            elif node.Id == 6:
                node_force_y = node.GetSolutionStepValue(DEM.CONTACT_FORCES_Y)
                expected_value = 150.1
                self.CheckForceY(node_force_y, expected_value, tolerance)

        super(DEM2D_ControlModuleTestSolution, self).Finalize()

    def CheckForceX(self, force, expected_value, tolerance):
        if abs(expected_value) > abs(force*tolerance) or abs(expected_value) < abs(force/tolerance):
            raise ValueError('Incorrect value for CONTACT_FORCES_X: expected value was '+ str(expected_value) + ' but received ' + str(force))

    def CheckForceY(self, force, expected_value, tolerance):
        if abs(expected_value) > abs(force*tolerance) or abs(expected_value) < abs(force/tolerance):
            raise ValueError('Incorrect value for CONTACT_FORCES_Y: expected value was '+ str(expected_value) + ' but received ' + str(force))


class TestDEM2DControlModule(KratosUnittest.TestCase):

    def setUp(self):
        pass

    def test_DEM2D_control_module(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM2D_control_module_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = KratosMultiphysics.Model()
        CreateAndRunStageInOneOpenMPThread(DEM2D_ControlModuleTestSolution, model, parameters_file_name)

    def tearDown(self):
        file_to_remove = os.path.join("DEM2D_control_module_tests_files", "TimesPartialRelease")
        kratos_utils.DeleteFileIfExisting(GetFilePath(file_to_remove))
        os.chdir(this_working_dir_backup)


if __name__ == "__main__":
    Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
