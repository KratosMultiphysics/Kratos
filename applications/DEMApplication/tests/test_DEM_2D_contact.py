import os
import KratosMultiphysics
from KratosMultiphysics import Logger
Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
import KratosMultiphysics.DEMApplication as DEM
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.DEMApplication.DEM_analysis_stage

import auxiliary_functions_for_tests

this_working_dir_backup = os.getcwd()

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class DEM2D_ContactTestSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage, KratosUnittest.TestCase):

    @classmethod
    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM2D_contact_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        tolerance = 1.001
        for node in self.rigid_face_model_part.Nodes:
            dem_pressure = node.GetSolutionStepValue(DEM.DEM_PRESSURE)
            contact_force = node.GetSolutionStepValue(DEM.CONTACT_FORCES_Y)
            if node.Id == 13:
                if self.time > 0.3:
                    expected_value = 45126
                    self.assertAlmostEqual(dem_pressure, expected_value, delta=tolerance)
                    expected_value = -23141
                    self.assertAlmostEqual(contact_force, expected_value, delta=tolerance)
            if node.Id == 22:
                if self.time > 0.3:
                    expected_value = 26712
                    self.assertAlmostEqual(dem_pressure, expected_value, delta=tolerance)
                    expected_value = -13698
                    self.assertAlmostEqual(contact_force, expected_value, delta=tolerance)


class TestDEM2DContact(KratosUnittest.TestCase):

    def setUp(self):
        pass

    @classmethod
    def test_DEM2D_contact(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM2D_contact_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = KratosMultiphysics.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(DEM2D_ContactTestSolution, model, parameters_file_name, 1)

    @staticmethod
    def tearDown():
        files_to_delete_list = []
        files_to_delete_list.append(os.path.join("TimesPartialRelease"))
        files_to_delete_list.append(os.path.join("DEM2D_contact_tests_files", "DEM2D_contact.post.lst"))
        files_to_delete_list.append(os.path.join("DEM2D_contact_tests_files", "flux_data_new.hdf5"))

        for to_erase_file in files_to_delete_list:
            if os.path.exists(to_erase_file):
                os.remove(to_erase_file)

        #............Getting rid of unuseful folders
        folders_to_delete_list      = []
        folders_to_delete_list.append(os.path.join("DEM2D_contact_tests_files", "DEM2D_contact_Graphs"))
        folders_to_delete_list.append(os.path.join("DEM2D_contact_tests_files", "DEM2D_contact_MPI_results"))
        folders_to_delete_list.append(os.path.join("DEM2D_contact_tests_files", "DEM2D_contact_Post_Files"))
        folders_to_delete_list.append(os.path.join("DEM2D_contact_tests_files", "DEM2D_contact_Results_and_Data"))


        for to_erase_folder in folders_to_delete_list:
            import shutil
            shutil.rmtree(to_erase_folder)

        os.chdir(this_working_dir_backup)

if __name__ == "__main__":
    Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
