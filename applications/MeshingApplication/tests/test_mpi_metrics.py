# We import the libraries
import KratosMultiphysics
import KratosMultiphysics.MeshingApplication as MeshingApplication

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.testing.utilities import ReadModelPart
import json
import os

def GetFilePath(fileName):
    return os.path.dirname(os.path.realpath(__file__)) + "/" + fileName

class TestMPIMetrics(KratosUnittest.TestCase):

    def test_mpi_hessian(self):
        KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

        # We create the model part
        current_model = KratosMultiphysics.Model()
        main_model_part = current_model.CreateModelPart("MainModelPart", 2)
        main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)

        # We add the variables needed
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)

        # We import the model main_model_part
        file_path = GetFilePath("/mmg_rectangle_test/remesh_rectangle")
        ReadModelPart(file_path, main_model_part)

        # We calculate the gradient of the distance variable
        find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(main_model_part)
        find_nodal_h.Execute()
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_AREA, 0.0, main_model_part.Nodes)

        ## Read reference
        file_name = "mmg_rectangle_test/rectangle_distance_values.json"
        distance_values_file_name = GetFilePath(file_name)

        with open(distance_values_file_name, 'r') as f:
            distance_values = json.load(f)

        for node in main_model_part.Nodes:
            distance_value = distance_values[str(node.Id)]
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, distance_value)

        metric_parameters = KratosMultiphysics.Parameters("""
        {
            "hessian_strategy_parameters"              :{
                "estimate_interpolation_error"     : false,
                "interpolation_error"              : 1.0e-6,
                "mesh_dependent_constant"          : 0.28125
            },
            "minimal_size"                     : 0.015,
            "maximal_size"                     : 0.5,
            "enforce_current"                      : false,
            "anisotropy_remeshing"                 : false,
            "enforce_anisotropy_relative_variable" : false
        }
        """)
        metric_process = MeshingApplication.ComputeHessianSolMetricProcess(main_model_part, KratosMultiphysics.DISTANCE, metric_parameters)
        metric_process.Execute()

        ## Read reference
        file_name = "mmg_rectangle_test/rectangle_pre_metric_result.json"
        reference_file_name = GetFilePath(file_name)
        with open(reference_file_name, 'r') as f:
            reference_values = json.load(f)

        for node in main_model_part.Nodes:
            metric_tensor = node.GetValue(MeshingApplication.METRIC_TENSOR_2D)
            ref_gradient = reference_values[str(node.Id)]
            for tensor_i, tensor_i_ref in zip(metric_tensor, ref_gradient):
                self.assertAlmostEqual(tensor_i, tensor_i_ref)

if __name__ == '__main__':
    KratosUnittest.main()
