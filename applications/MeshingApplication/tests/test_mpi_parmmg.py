import KratosMultiphysics
from KratosMultiphysics import ParallelEnvironment, IsDistributedRun
import KratosMultiphysics.MeshingApplication
from KratosMultiphysics.compare_two_files_check_process import CompareTwoFilesCheckProcess
import KratosMultiphysics.kratos_utilities as kratos_utilities
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.testing.utilities import ReadModelPart
import os
import math

def GetFilePath(fileName):
    return os.path.dirname(os.path.realpath(__file__)) + "/" + fileName

class TestMPIParMmg(KratosUnittest.TestCase):

    @KratosUnittest.skipUnless(IsDistributedRun() and ParallelEnvironment.GetDefaultSize() == 2, "Test designed to be run with two ranks.")
    def test_mpi_sphere(self):
        KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

        # We create the model part
        current_model = KratosMultiphysics.Model()
        main_model_part = current_model.CreateModelPart("MainModelPart")
        main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)

        # We add the variables needed
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)


        # We import the model main_model_part
        file_path = GetFilePath("/parmmg_eulerian_test/background_mesh_sphere")
        ReadModelPart(file_path, main_model_part)

        communicator = main_model_part.GetCommunicator().GetDataCommunicator()

        for node in main_model_part.Nodes:
            distance = math.sqrt(node.X**2+node.Y**2+node.Z**2) - 1.0/2.0
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE,distance)

        ##COMPUTE DISTANCE GRADIENT AND NODAL_H
        local_gradient = KratosMultiphysics.ComputeNodalGradientProcess3D(main_model_part,
        KratosMultiphysics.DISTANCE,
        KratosMultiphysics.DISTANCE_GRADIENT,
        KratosMultiphysics.NODAL_AREA)

        local_gradient.Execute()
        find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(main_model_part)
        find_nodal_h.Execute()

        ##COMPUTE LEVEL SET METRIC
        metric_parameters = KratosMultiphysics.Parameters("""
        {
            "minimal_size"                             : 1.5,
            "sizing_parameters": {
                "reference_variable_name"               : "DISTANCE",
                "boundary_layer_max_distance"           : 2.0,
                "interpolation"                         : "constant"
            },
            "enforce_current"                      : false,
            "anisotropy_remeshing"                 : false
        }
        """)
        metric_process = KratosMultiphysics.MeshingApplication.ComputeLevelSetSolMetricProcess3D(
            main_model_part,
            KratosMultiphysics.DISTANCE_GRADIENT,
            metric_parameters)
        metric_process.Execute()

        ##PERFORM REMESHING
        pmmg_parameters = KratosMultiphysics.Parameters("""
        {
            "filename"                         : "output",
            "save_external_files"              : true,
            "save_colors_files"                : true,
            "initialize_entities"              : false,
            "preserve_flags"                   : false
        }
        """)
        pmmg_parameters["filename"].SetString(GetFilePath(pmmg_parameters["filename"].GetString()))
        pmmg_process = KratosMultiphysics.MeshingApplication.ParMmgProcess3D(main_model_part.GetRootModelPart(), pmmg_parameters)
        pmmg_process.Execute()

        # Checking color element and condition maps.
        check_parameters = KratosMultiphysics.Parameters("""
        {
            "reference_file_name"   : "parmmg_eulerian_test/cond_ref_map.json",
            "output_file_name"      : "output_file_name",
            "comparison_type"       : "deterministic"

        }
        """)

        check_parameters["output_file_name"].SetString(GetFilePath("output_step=0_"+str(communicator.Rank())+".cond.ref.json"))
        check_parameters["reference_file_name"].SetString(GetFilePath(check_parameters["reference_file_name"].GetString()))
        check_files = CompareTwoFilesCheckProcess(check_parameters)
        check_files.ExecuteInitialize()
        check_files.ExecuteBeforeSolutionLoop()
        check_files.ExecuteInitializeSolutionStep()
        check_files.ExecuteFinalizeSolutionStep()
        check_files.ExecuteFinalize()

        check_parameters = KratosMultiphysics.Parameters("""
        {
            "reference_file_name"   : "parmmg_eulerian_test/elem_ref_map.json",
            "output_file_name"      : "output_file_name",
            "comparison_type"       : "deterministic"

        }
        """)

        check_parameters["output_file_name"].SetString(GetFilePath("output_step=0_"+str(communicator.Rank())+".elem.ref.json"))
        check_parameters["reference_file_name"].SetString(GetFilePath(check_parameters["reference_file_name"].GetString()))
        check_files = CompareTwoFilesCheckProcess(check_parameters)
        check_files.ExecuteInitialize()
        check_files.ExecuteBeforeSolutionLoop()
        check_files.ExecuteInitializeSolutionStep()
        check_files.ExecuteFinalizeSolutionStep()
        check_files.ExecuteFinalize()

        pmmg_process.OutputMdpa()

        # We check the solution
        check_parameters = KratosMultiphysics.Parameters("""
        {
            "reference_file_name"   : "reference_file_name",
            "output_file_name"      : "output_file_name",
            "comparison_type"       : "deterministic"

        }
        """)
        check_parameters["reference_file_name"].SetString(GetFilePath("parmmg_eulerian_test/parmmg_sphere_reference_mdpa_rank_"+str(communicator.Rank())+".mdpa"))
        check_parameters["output_file_name"].SetString(GetFilePath("output_"+str(communicator.Rank())+".mdpa"))
        check_files = CompareTwoFilesCheckProcess(check_parameters)

        check_files.ExecuteInitialize()
        check_files.ExecuteBeforeSolutionLoop()
        check_files.ExecuteInitializeSolutionStep()
        check_files.ExecuteFinalizeSolutionStep()
        check_files.ExecuteFinalize()

        for file_name in os.listdir(GetFilePath("")):
            if file_name.endswith(".json") or file_name.endswith(".mdpa") or file_name.endswith(".mesh") or  file_name.endswith(".sol"):
                kratos_utilities.DeleteFileIfExisting(GetFilePath(file_name))
        kratos_utilities.DeleteTimeFiles(os.getcwd())



if __name__ == '__main__':
    KratosUnittest.main()
