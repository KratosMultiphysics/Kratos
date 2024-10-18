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

class TestClusters(KratosUnittest.TestCase):

    def setUp(self):
        source = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "custom_elements", "custom_clusters", "linecluster3D.clu")
        destination = os.path.join(os.path.dirname(os.path.realpath(__file__)), subfolder_name, "linecluster3D.clu")
        copyfile(source, destination)

    def test_clusters_1(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), subfolder_name)
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        with open(parameters_file_name,'r') as parameter_file:
            project_parameters = KratosMultiphysics.Parameters(parameter_file.read())
        model = KratosMultiphysics.Model()

        with auxiliary_functions_for_tests.controlledExecutionScope(path):
            modified_dem_analysis_stage = KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage(model, project_parameters)

            backup_original_finalize = modified_dem_analysis_stage.Finalize
            def ModifiedFinalize(*args, **kwargs):
                def Aux(fnc, *args, **kwargs):
                    tol = 1.0e-3
                    node = modified_dem_analysis_stage.spheres_model_part.GetNode(22)
                    velz = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z)
                    self.assertAlmostEqual(velz, -0.3888113323025093, delta=tol)

                    node = modified_dem_analysis_stage.spheres_model_part.GetNode(31)
                    velz = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z)
                    self.assertAlmostEqual(velz, 0.029044273544728355, delta=tol)

                    fnc(*args, **kwargs)

                Aux(backup_original_finalize, *args, **kwargs)

            modified_dem_analysis_stage.Finalize = ModifiedFinalize

            auxiliary_functions_for_tests.RunStageInSelectedNumberOfOpenMPThreads(modified_dem_analysis_stage, 1)

    def test_clusters_2(self):
        self.check_mark_1 = False
        self.check_mark_2 = False
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), subfolder_name)
        parameters_file_name = os.path.join(path, "ProjectParametersDEM2.json")
        with open(parameters_file_name,'r') as parameter_file:
            project_parameters = KratosMultiphysics.Parameters(parameter_file.read())
        model = KratosMultiphysics.Model()

        with auxiliary_functions_for_tests.controlledExecutionScope(path):
            modified_dem_analysis_stage = KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage(model, project_parameters)
            backup_original_finalize_solution_step = modified_dem_analysis_stage.FinalizeSolutionStep

            def ModifiedFinalizeSolutionStep(*args, **kwargs):
                def Aux(fnc, *args, **kwargs):
                    self.assertTrue(modified_dem_analysis_stage.DEM_parameters["BoundingBoxOption"].GetBool())
                    self.assertAlmostEqual(modified_dem_analysis_stage.bounding_box_time_limits[0], 1.0e-6, delta=1e-12)
                    self.assertAlmostEqual(modified_dem_analysis_stage.bounding_box_time_limits[1], 1000.0, delta=1e-12)
                    if modified_dem_analysis_stage.time < 4.4e-6:
                        self.assertEqual(modified_dem_analysis_stage.spheres_model_part.NumberOfNodes(0), 128)
                        self.assertEqual(modified_dem_analysis_stage.spheres_model_part.NumberOfElements(0), 128)
                        self.assertEqual(modified_dem_analysis_stage.cluster_model_part.NumberOfNodes(0), 2)
                        self.assertEqual(modified_dem_analysis_stage.cluster_model_part.NumberOfElements(0), 2)
                        self.check_mark_1 = True
                    else:
                        self.assertEqual(modified_dem_analysis_stage.spheres_model_part.NumberOfNodes(0), 64)
                        self.assertEqual(modified_dem_analysis_stage.spheres_model_part.NumberOfElements(0), 64)
                        self.assertEqual(modified_dem_analysis_stage.cluster_model_part.NumberOfNodes(0), 1)
                        self.assertEqual(modified_dem_analysis_stage.cluster_model_part.NumberOfElements(0), 1)
                        self.check_mark_2 = True
                    fnc(*args, **kwargs)
                Aux(backup_original_finalize_solution_step, *args, **kwargs)

            backup_original_finalize = modified_dem_analysis_stage.Finalize
            def ModifiedFinalize(*args, **kwargs):
                def Aux(fnc, *args, **kwargs):
                    self.assertTrue(self.check_mark_1)
                    self.assertTrue(self.check_mark_2)
                    fnc(*args, **kwargs)
                Aux(backup_original_finalize, *args, **kwargs)


            modified_dem_analysis_stage.FinalizeSolutionStep = ModifiedFinalizeSolutionStep
            modified_dem_analysis_stage.Finalize = ModifiedFinalize

            auxiliary_functions_for_tests.RunStageInSelectedNumberOfOpenMPThreads(modified_dem_analysis_stage, 1)

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
