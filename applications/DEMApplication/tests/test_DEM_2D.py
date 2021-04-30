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

class DEM2DTestSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage, KratosUnittest.TestCase):

    @classmethod
    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM2D_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        tolerance = 1e-3
        for node in self.spheres_model_part.Nodes:
            normal_impact_vel = node.GetSolutionStepValue(Kratos.VELOCITY_X)
            if node.Id == 1:
                if self.time > 0.2:
                    self.assertAlmostEqual(normal_impact_vel, 6.076801447242313, delta=tolerance)
            if node.Id == 2:
                if self.time > 0.2:
                    self.assertAlmostEqual(normal_impact_vel, 8.604163136887411, delta=tolerance)
            if node.Id == 3:
                if self.time > 0.2:
                    self.assertAlmostEqual(normal_impact_vel, 10.016439272775422, delta=tolerance)

class TestDEM2D(KratosUnittest.TestCase):

    def setUp(self):
        pass

    @classmethod
    def test_DEM2D_1(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM2D_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = Kratos.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(DEM2DTestSolution, model, parameters_file_name, 1)


    @staticmethod
    def tearDown():
        files_to_delete_list = []
        files_to_delete_list.append(os.path.join("TimesPartialRelease"))
        files_to_delete_list.append(os.path.join("DEM2D_tests_files", "impact2D.post.lst"))
        files_to_delete_list.append(os.path.join("DEM2D_tests_files", "flux_data_new.hdf5"))

        for to_erase_file in files_to_delete_list:
            if os.path.exists(to_erase_file):
                os.remove(to_erase_file)

        #............Getting rid of unuseful folders
        folders_to_delete_list      = []
        folders_to_delete_list.append(os.path.join("DEM2D_tests_files", "impact2D_Graphs"))
        folders_to_delete_list.append(os.path.join("DEM2D_tests_files", "impact2D_MPI_results"))
        folders_to_delete_list.append(os.path.join("DEM2D_tests_files", "impact2D_Post_Files"))
        folders_to_delete_list.append(os.path.join("DEM2D_tests_files", "impact2D_Results_and_Data"))


        for to_erase_folder in folders_to_delete_list:
            import shutil
            shutil.rmtree(to_erase_folder)

        os.chdir(this_working_dir_backup)


if __name__ == "__main__":
    Kratos.Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
