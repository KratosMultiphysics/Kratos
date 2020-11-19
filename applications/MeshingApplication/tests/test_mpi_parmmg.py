import KratosMultiphysics
from KratosMultiphysics import ParallelEnvironment, IsDistributedRun
import KratosMultiphysics.MeshingApplication
import KratosMultiphysics.kratos_utilities as kratos_utilities
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.testing.utilities import ReadModelPart
from KratosMultiphysics.testing.utilities import ReadSerialModelPart
import os
import math
import json

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
            "preserve_flags"                   : false,
            "echo_level"                       : 0
        }
        """)
        pmmg_parameters["filename"].SetString(GetFilePath(pmmg_parameters["filename"].GetString()))
        pmmg_process = KratosMultiphysics.MeshingApplication.ParMmgProcess3D(main_model_part.GetRootModelPart(), pmmg_parameters)
        pmmg_process.Execute()

        reference_file_name = GetFilePath("parmmg_eulerian_test/cond_ref_map.json")
        result_file_name = GetFilePath("output_step=0_"+str(communicator.Rank())+".cond.ref.json")
        self._CompareColorFiles(reference_file_name, result_file_name)

        reference_file_name = GetFilePath("parmmg_eulerian_test/elem_ref_map.json")
        result_file_name = GetFilePath("output_step=0_"+str(communicator.Rank())+".elem.ref.json")
        self._CompareColorFiles(reference_file_name, result_file_name)

        ref_mdpa_filename = GetFilePath("parmmg_eulerian_test/parmmg_sphere_reference_mdpa_rank_"+str(communicator.Rank()))
        ref_model_part = current_model.CreateModelPart("Reference")
        # We add the variables needed
        ref_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        ref_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)
        ref_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)
        ReadSerialModelPart(ref_mdpa_filename, ref_model_part)

        self._CheckModelPart(ref_model_part, main_model_part)
        for ref_sub_model_part in ref_model_part.SubModelParts:
            sub_model_part_name = ref_sub_model_part.Name
            if main_model_part.HasSubModelPart(sub_model_part_name):
                result_sub_model_part = main_model_part.GetSubModelPart(sub_model_part_name)
            else:
                raise(Exception("Submodelpart", sub_model_part_name, "does not exist in the remeshed model part"))
            self._CheckModelPart(ref_sub_model_part, result_sub_model_part)

        for file_name in os.listdir(GetFilePath("")):
            if file_name.endswith(".json") or file_name.endswith(".mdpa") or file_name.endswith(".mesh") or  file_name.endswith(".sol"):
                kratos_utilities.DeleteFileIfExisting(GetFilePath(file_name))
        kratos_utilities.DeleteTimeFiles(os.getcwd())

    def _CompareColorFiles(self, ref_dict_filename, result_dict_file_name):

        with open(ref_dict_filename, 'r') as f:
            reference_values = json.load(f)

        with open(result_dict_file_name, 'r') as f:
            result_values = json.load(f)

        self.assertEqual(len(reference_values.keys()), len(result_values.keys()))
        for key_ref, key_result in zip(reference_values.keys(), result_values.keys()):
            self.assertEqual(reference_values[key_ref], result_values[key_result])

    def _CheckModelPart(self, ref_model_part, result_model_part):
        self.assertEqual(ref_model_part.NumberOfNodes(), result_model_part.NumberOfNodes())
        self.assertEqual(ref_model_part.NumberOfElements(), result_model_part.NumberOfElements())
        self.assertEqual(ref_model_part.NumberOfConditions(), result_model_part.NumberOfConditions())

if __name__ == '__main__':
    KratosUnittest.main()
