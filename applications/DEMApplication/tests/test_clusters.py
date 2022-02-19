import os
import KratosMultiphysics
from KratosMultiphysics import Logger
Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
import KratosMultiphysics.DEMApplication as DEM
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.DEMApplication.DEM_analysis_stage

import KratosMultiphysics.kratos_utilities as kratos_utils

import auxiliary_functions_for_tests

from shutil import copyfile

this_working_dir_backup = os.getcwd()

subfolder_name = "cluster_tests_files"

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class ClustersTestSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage, KratosUnittest.TestCase):

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def Finalize(self):
        tol = 1e-3
        for node in self.spheres_model_part.Nodes:
            if node.Id == 21:
                velz = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z)
                self.assertAlmostEqual(velz, -0.3888113323025093, delta=tol)
            if node.Id == 30:
                velz = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z)
                self.assertAlmostEqual(velz, 0.029044273544728355, delta=tol)

        del node

        super().Finalize()


class TestClusters(KratosUnittest.TestCase):

    def setUp(self):
        source = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "custom_elements", "custom_clusters", "linecluster3D.clu")
        destination = os.path.join(os.path.dirname(os.path.realpath(__file__)), subfolder_name, "linecluster3D.clu")
        copyfile(source, destination)

    @classmethod
    def test_clusters_1(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), subfolder_name)
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = KratosMultiphysics.Model()
        with auxiliary_functions_for_tests.controlledExecutionScope(path):
            auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(ClustersTestSolution, model, parameters_file_name, 1)

    def tearDown(self):
        file_to_remove = os.path.join(subfolder_name, "TimesPartialRelease")
        kratos_utils.DeleteFileIfExisting(GetFilePath(file_to_remove))
        file_to_remove = os.path.join(subfolder_name, "flux_data_new.hdf5")
        kratos_utils.DeleteFileIfExisting(GetFilePath(file_to_remove))
        file_to_remove = os.path.join(subfolder_name, "linecluster3D.clu")
        kratos_utils.DeleteFileIfExisting(GetFilePath(file_to_remove))
        os.chdir(this_working_dir_backup)


if __name__ == "__main__":
    Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
