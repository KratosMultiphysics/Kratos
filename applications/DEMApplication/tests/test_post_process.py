import os
import shutil
import KratosMultiphysics
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

class TestPostProcessClass1(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage, KratosUnittest.TestCase):

    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "post_process_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def Finalize(self):
        self.procedures.RemoveFoldersWithResults(str(self.main_path), str(self.problem_name), '')
        super().Finalize()

class TestPostProcess(KratosUnittest.TestCase):

    def setUp(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "post_process_tests_files")
        tests_path = os.path.dirname(os.path.realpath(__file__))
        clusters_path = os.path.join(tests_path, "..", "custom_elements", "custom_clusters")
        ring_cluster_path = os.path.join(clusters_path, "ringcluster3D.clu")
        shutil.copy(ring_cluster_path, path)

    def test_gid_printing_many_results(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "post_process_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = KratosMultiphysics.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(TestPostProcessClass1, model, parameters_file_name, auxiliary_functions_for_tests.GetHardcodedNumberOfThreads())

    def tearDown(self):
        file_to_remove = os.path.join("post_process_tests_files", "ringcluster3D.clu")
        kratos_utils.DeleteFileIfExisting(GetFilePath(file_to_remove))
        os.chdir(this_working_dir_backup)


if __name__ == "__main__":
    Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
