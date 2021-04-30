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

class DEM2D_InletTestSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage, KratosUnittest.TestCase):

    @classmethod
    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM2D_inlet_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        tolerance = 1.001
        for node in self.spheres_model_part.Nodes:
            node_vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)
            node_force = node.GetSolutionStepValue(KratosMultiphysics.TOTAL_FORCES_Y)
            if node.Id == 6:
                if self.time >= 1.15:
                    Logger.PrintInfo(node_vel)
                    Logger.PrintInfo(node_force)
                    self.assertAlmostEqual(node_vel, 0.380489240, delta=tolerance)
                    self.assertAlmostEqual(node_force, -120983.1002, delta=tolerance)

class TestDEM2DInlet(KratosUnittest.TestCase):

    def setUp(self):
        pass

    @classmethod
    def test_DEM2D_inlet(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM2D_inlet_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = KratosMultiphysics.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(DEM2D_InletTestSolution, model, parameters_file_name, 1)

    def tearDown(self):
        file_to_remove = os.path.join("DEM2D_inlet_tests_files", "TimesPartialRelease")
        kratos_utils.DeleteFileIfExisting(GetFilePath(file_to_remove))
        file_to_remove = os.path.join("DEM2D_inlet_tests_files", "flux_data_new.hdf5")
        kratos_utils.DeleteFileIfExisting(GetFilePath(file_to_remove))
        os.chdir(this_working_dir_backup)


    @staticmethod
    def tearDown():
        files_to_delete_list = []
        files_to_delete_list.append(os.path.join("TimesPartialRelease"))
        files_to_delete_list.append(os.path.join("DEM2D_inlet_tests_files", "DEM2D_inlet.post.lst"))
        files_to_delete_list.append(os.path.join("DEM2D_inlet_tests_files", "flux_data_new.hdf5"))

        for to_erase_file in files_to_delete_list:
            if os.path.exists(to_erase_file):
                os.remove(to_erase_file)

        #............Getting rid of unuseful folders
        folders_to_delete_list      = []
        folders_to_delete_list.append(os.path.join("DEM2D_inlet_tests_files", "DEM2D_inlet_Graphs"))
        folders_to_delete_list.append(os.path.join("DEM2D_inlet_tests_files", "DEM2D_inlet_MPI_results"))
        folders_to_delete_list.append(os.path.join("DEM2D_inlet_tests_files", "DEM2D_inlet_Post_Files"))
        folders_to_delete_list.append(os.path.join("DEM2D_inlet_tests_files", "DEM2D_inlet_Results_and_Data"))


        for to_erase_folder in folders_to_delete_list:
            import shutil
            shutil.rmtree(to_erase_folder)

        os.chdir(this_working_dir_backup)



if __name__ == "__main__":
    Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
