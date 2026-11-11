import os
import KratosMultiphysics as Kratos
import KratosMultiphysics.DEMApplication as DEM
from KratosMultiphysics import Logger
Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.DEMApplication.DEM_analysis_stage
import KratosMultiphysics.kratos_utilities as kratos_utils

import auxiliary_functions_for_tests

this_working_dir_backup = os.getcwd()

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class SmoothJointModelTestSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage, KratosUnittest.TestCase):

    @classmethod
    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "dem_3d_smooth_joint_model_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()
        self.spheres_model_part.ProcessInfo[DEM.IS_TIME_TO_UPDATE_CONTACT_ELEMENT] = True

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        tolerance = 1e-7
        for node in self.spheres_model_part.Nodes:
            velocity = node.GetSolutionStepValue(Kratos.VELOCITY)
            angular_velocity = node.GetSolutionStepValue(Kratos.ANGULAR_VELOCITY)
            if node.Id == 1:
                if self.time > 0.00098 and self.time < 0.000982:
                    expected_value = 0.0
                    self.CheckValues(velocity, 0, expected_value, tolerance)
                    expected_value = 0.0
                    self.CheckValues(velocity, 1, expected_value, tolerance)
                    expected_value = 0.0
                    self.CheckValues(velocity, 2, expected_value, tolerance)
                    expected_value = 0.0
                    self.CheckValues(angular_velocity, 0, expected_value, tolerance)
                    expected_value = 0.0
                    self.CheckValues(angular_velocity, 1, expected_value, tolerance)
                    expected_value = 0.0
                    self.CheckValues(angular_velocity, 2, expected_value, tolerance)

    def CheckValues(self, variable, component, expected_value, tolerance):
        self.assertAlmostEqual(variable[component], expected_value, delta=tolerance)

    def Finalize(self):
        self.procedures.RemoveFoldersWithResults(str(self.main_path), str(self.problem_name), '')
        super().Finalize()

class SmoothJointMirroredContactsTestSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage, KratosUnittest.TestCase):

    def __init__(self, model, parameters):
        super().__init__(model, parameters)
        self.mirrored_contacts_checked = False

    @classmethod
    def GetMainPath(cls):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "dem_3d_smooth_joint_mirrored_contacts_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

        if self.time <= 5.0 * self.DEM_parameters["MaxTimeStep"].GetDouble():
            return

        contact_elements = list(self.contact_model_part.Elements)
        if len(contact_elements) != 2:
            return

        self._GetSolver().PrepareContactElementsForPrinting()

        contacts_by_nodes = {
            tuple(sorted(node.Id for node in element.GetGeometry())): element
            for element in contact_elements
        }
        self.assertTrue(set(contacts_by_nodes) == {(1, 2), (3, 4)})

        contact_sigmas = []
        contact_taus = []
        for node_pair in ((1, 2), (3, 4)):
            contact = contacts_by_nodes[node_pair]
            contact_sigma = contact.GetValue(DEM.CONTACT_SIGMA)
            contact_tau = contact.GetValue(DEM.CONTACT_TAU)
            contact_failure = contact.GetValue(DEM.CONTACT_FAILURE)

            self.assertGreater(contact_sigma, 0.0)
            self.assertGreater(contact_tau, 0.0)
            self.assertTrue(contact_failure == 0.0)
            contact_sigmas.append(contact_sigma)
            contact_taus.append(contact_tau)

        # The two contacts have the same normal closure but their particle
        # radii (and hence their joint areas) differ by a factor of four.
        # Stiffness properties are specified per unit area [Pa/m], so both
        # contacts must develop the same stress. This also guards against
        # accidentally treating JOINT_NORMAL_STIFFNESS as a contact stiffness
        # in N/m, which would make stress inversely proportional to area.
        self.assertAlmostEqual(contact_sigmas[0], contact_sigmas[1], delta=1.0e-6)
        self.assertAlmostEqual(contact_taus[0], contact_taus[1], delta=1.0e-6)
        self.mirrored_contacts_checked = True

    def Finalize(self):
        self.assertTrue(self.mirrored_contacts_checked)
        super().Finalize()

class SmoothJointDilationTestSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage, KratosUnittest.TestCase):

    def __init__(self, model, parameters):
        super().__init__(model, parameters)
        self.dilation_checked = False

    @classmethod
    def GetMainPath(cls):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "dem_3d_smooth_joint_dilation_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

        if self.time <= 8.0 * self.DEM_parameters["MaxTimeStep"].GetDouble():
            return

        contact_elements = list(self.contact_model_part.Elements)
        if len(contact_elements) != 2:
            return

        self._GetSolver().PrepareContactElementsForPrinting()
        contacts_by_nodes = {
            tuple(sorted(node.Id for node in element.GetGeometry())): element
            for element in contact_elements
        }
        self.assertTrue(set(contacts_by_nodes) == {(1, 2), (3, 4)})

        contact_without_dilation = contacts_by_nodes[(1, 2)]
        contact_with_dilation = contacts_by_nodes[(3, 4)]

        sigma_without_dilation = contact_without_dilation.GetValue(DEM.CONTACT_SIGMA)
        sigma_with_dilation = contact_with_dilation.GetValue(DEM.CONTACT_SIGMA)
        tau_without_dilation = contact_without_dilation.GetValue(DEM.CONTACT_TAU)
        tau_with_dilation = contact_with_dilation.GetValue(DEM.CONTACT_TAU)

        self.assertTrue(contact_without_dilation.GetValue(DEM.CONTACT_FAILURE) == 2.0)
        self.assertTrue(contact_with_dilation.GetValue(DEM.CONTACT_FAILURE) == 2.0)
        self.assertAlmostEqual(tau_without_dilation, 0.0, delta=1.0e-8)
        self.assertAlmostEqual(tau_with_dilation, 0.0, delta=1.0e-8)
        self.assertGreater(sigma_without_dilation, 0.0)
        self.assertGreater(sigma_with_dilation, sigma_without_dilation)
        self.dilation_checked = True

    def Finalize(self):
        self.assertTrue(self.dilation_checked)
        super().Finalize()

class TestSmoothJointModel(KratosUnittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "dem_3d_smooth_joint_mirrored_contacts_tests_files")
        kratos_utils.DeleteDirectoryIfExisting(os.path.join(path, "test_SJM_mirrored_Post_Files"))
        kratos_utils.DeleteDirectoryIfExisting(os.path.join(path, "test_SJM_mirrored_Graphs"))
        kratos_utils.DeleteFileIfExisting(os.path.join(path, "test_SJM_mirrored.post.lst"))
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "dem_3d_smooth_joint_dilation_tests_files")
        kratos_utils.DeleteDirectoryIfExisting(os.path.join(path, "test_SJM_dilation_Post_Files"))
        kratos_utils.DeleteDirectoryIfExisting(os.path.join(path, "test_SJM_dilation_Graphs"))
        kratos_utils.DeleteFileIfExisting(os.path.join(path, "test_SJM_dilation.post.lst"))

    @classmethod
    def test_SmoothJointModel_1(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "dem_3d_smooth_joint_model_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = Kratos.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(SmoothJointModelTestSolution, model, parameters_file_name, 1)

    def test_SmoothJointModel_MirroredContacts(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "dem_3d_smooth_joint_mirrored_contacts_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = Kratos.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(SmoothJointMirroredContactsTestSolution, model, parameters_file_name, 1)

    def test_SmoothJointModel_Dilation(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "dem_3d_smooth_joint_dilation_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = Kratos.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(SmoothJointDilationTestSolution, model, parameters_file_name, 1)

if __name__ == "__main__":
    Kratos.Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
