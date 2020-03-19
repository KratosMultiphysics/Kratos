import os
import KratosMultiphysics
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
    omp_utils = KratosMultiphysics.OpenMPUtils()
    if "OMP_NUM_THREADS" in os.environ:
        initial_number_of_threads = os.environ['OMP_NUM_THREADS']
        omp_utils.SetNumThreads(1)

    with open(parameters_file_name,'r') as parameter_file:
        project_parameters = KratosMultiphysics.Parameters(parameter_file.read())

    my_obj(model, project_parameters).Run()

    if "OMP_NUM_THREADS" in os.environ:
        omp_utils.SetNumThreads(int(initial_number_of_threads))

class DEM3D_ContactTestSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage):

    @classmethod
    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM3D_contact_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def FinalizeTimeStep(self, time):
        tolerance = 1.001
        for node in self.rigid_face_model_part.Nodes:
            dem_pressure = node.GetSolutionStepValue(DEM.DEM_PRESSURE)
            contact_force = node.GetSolutionStepValue(DEM.CONTACT_FORCES_Z)
            if node.Id == 9:
                if time > 0.35:
                    expected_value = 1621
                    self.CheckPressure(dem_pressure, expected_value, tolerance)
                    expected_value = -6484
                    self.CheckContactF(contact_force, expected_value, tolerance)
            if node.Id == 13:
                if time > 0.35:
                    expected_value = 841
                    self.CheckPressure(dem_pressure, expected_value, tolerance)
                    expected_value = -3366
                    self.CheckContactF(contact_force, expected_value, tolerance)



    @classmethod
    def CheckPressure(self, dem_pressure, expected_value, tolerance):
        if expected_value > dem_pressure*tolerance or expected_value < dem_pressure/tolerance:
            raise ValueError('Incorrect value for DEM_PRESSURE: expected value was '+ str(expected_value) + ' but received ' + str(dem_pressure))

    @classmethod
    def CheckContactF(self, contact_force, expected_value, tolerance):
        if abs(expected_value) > abs(contact_force*tolerance) or abs(expected_value) < abs(contact_force/tolerance):
            raise ValueError('Incorrect value for CONTACT_FORCES_X: expected value was '+ str(expected_value) + ' but received ' + str(contact_force))

    def Finalize(self):
        super(DEM3D_ContactTestSolution, self).Finalize()



class TestDEM3DContact(KratosUnittest.TestCase):

    def setUp(self):
        pass

    @classmethod
    def test_DEM3D_contact(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM3D_contact_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = KratosMultiphysics.Model()
        CreateAndRunStageInOneOpenMPThread(DEM3D_ContactTestSolution, model, parameters_file_name)

    def tearDown(self):
        file_to_remove = os.path.join("DEM3D_contact_tests_files", "TimesPartialRelease")
        kratos_utils.DeleteFileIfExisting(GetFilePath(file_to_remove))
        os.chdir(this_working_dir_backup)


if __name__ == "__main__":
    Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
