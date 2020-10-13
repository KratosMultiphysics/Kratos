import os
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics
from KratosMultiphysics.testing.utilities import ReadModelPart
import json

def GetFilePath(fileName):
    return os.path.dirname(os.path.realpath(__file__)) + "/" + fileName

class TestMPIProcesses(KratosUnittest.TestCase):

    def testComputeNodalGradientProcess(self):

        current_model = KratosMultiphysics.Model()
        main_model_part = current_model.CreateModelPart("main_model_part")
        main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]=2

        ## Add variables to the model part
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)

        ## Serial partition of the original .mdpa file
        input_filename = GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_mpi_serializer")
        ReadModelPart(input_filename, main_model_part)

        ParallelFillCommunicator = KratosMultiphysics.mpi.ParallelFillCommunicator(main_model_part)
        ParallelFillCommunicator.Execute()

        for node in main_model_part.Nodes:
            distance = node.X**2+node.Y**2+node.Z**2 - 1
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE,distance)
            node.SetValue(KratosMultiphysics.DISTANCE,distance)

        ## Read reference
        file_name = "auxiliar_files_for_python_unittest/reference_files/test_compute_nodal_gradient_process_results.json"
        reference_file_name = GetFilePath(file_name)
        with open(reference_file_name, 'r') as f:
            reference_values = json.load(f)

        gradient_process_hist_hist = KratosMultiphysics.ComputeNodalGradientProcess(main_model_part,
        KratosMultiphysics.DISTANCE,
        KratosMultiphysics.DISTANCE_GRADIENT,
        KratosMultiphysics.NODAL_AREA)
        gradient_process_hist_hist.Execute()

        for node in main_model_part.Nodes:
                distance_gradient = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE_GRADIENT)
                ref_gradient = reference_values[str(node.Id)]
                for gradient_i, gradient_i_ref in zip(distance_gradient, ref_gradient):
                    self.assertAlmostEqual(gradient_i, gradient_i_ref)

        non_historical_origin_variable = True
        gradient_process_non_hist_hist = KratosMultiphysics.ComputeNodalGradientProcess(main_model_part,
        KratosMultiphysics.DISTANCE,
        KratosMultiphysics.DISTANCE_GRADIENT,
        KratosMultiphysics.NODAL_AREA,
        non_historical_origin_variable)
        gradient_process_non_hist_hist.Execute()

        for node in main_model_part.Nodes:
            distance_gradient = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE_GRADIENT)
            ref_gradient = reference_values[str(node.Id)]
            for gradient_i, gradient_i_ref in zip(distance_gradient, ref_gradient):
                self.assertAlmostEqual(gradient_i, gradient_i_ref)

        gradient_parameters = KratosMultiphysics.Parameters(""" {
            "origin_variable"                : "DISTANCE",
            "gradient_variable"              : "DISTANCE_GRADIENT",
            "non_historical_origin_variable" :  false
        }
        """)
        gradient_process = KratosMultiphysics.ComputeNodalGradientProcess(main_model_part,
        gradient_parameters)
        gradient_process.Execute()

        for node in main_model_part.Nodes:
            distance_gradient = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE_GRADIENT)
            ref_gradient = reference_values[str(node.Id)]
            for i_gradient, i_reference in zip(distance_gradient, ref_gradient):
                self.assertAlmostEqual(i_gradient, i_reference)

        gradient_parameters = KratosMultiphysics.Parameters(""" {
            "origin_variable"                : "DISTANCE",
            "gradient_variable"              : "DISTANCE_GRADIENT",
            "non_historical_origin_variable" :  true
        }
        """)
        gradient_process = KratosMultiphysics.ComputeNodalGradientProcess(main_model_part,
        gradient_parameters)
        gradient_process.Execute()

        for node in main_model_part.Nodes:
            distance_gradient = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE_GRADIENT)
            ref_gradient = reference_values[str(node.Id)]
            for i_gradient, i_reference in zip(distance_gradient, ref_gradient):
                self.assertAlmostEqual(i_gradient, i_reference)

    def testComputeNonHistoricalNodalGradientProcess(self):

        current_model = KratosMultiphysics.Model()
        main_model_part = current_model.CreateModelPart("main_model_part")
        main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]=2

        ## Add variables to the model part
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)

        ## Serial partition of the original .mdpa file
        input_filename = GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_mpi_serializer")
        ReadModelPart(input_filename, main_model_part)

        ParallelFillCommunicator = KratosMultiphysics.mpi.ParallelFillCommunicator(main_model_part)
        ParallelFillCommunicator.Execute()

        for node in main_model_part.Nodes:
            distance = node.X**2+node.Y**2+node.Z**2 - 1
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE,distance)
            node.SetValue(KratosMultiphysics.DISTANCE,distance)

        ## Read reference
        file_name = "auxiliar_files_for_python_unittest/reference_files/test_compute_nodal_gradient_process_results.json"
        reference_file_name = GetFilePath(file_name)
        with open(reference_file_name, 'r') as f:
            reference_values = json.load(f)

        gradient_process_hist_non_hist = KratosMultiphysics.ComputeNonHistoricalNodalGradientProcess(main_model_part,
        KratosMultiphysics.DISTANCE,
        KratosMultiphysics.DISTANCE_GRADIENT,
        KratosMultiphysics.NODAL_AREA)
        gradient_process_hist_non_hist.Execute()

        for node in main_model_part.Nodes:
            distance_gradient = node.GetValue(KratosMultiphysics.DISTANCE_GRADIENT)
            ref_gradient = reference_values[str(node.Id)]
            for gradient_i, gradient_i_ref in zip(distance_gradient, ref_gradient):
                self.assertAlmostEqual(gradient_i, gradient_i_ref)

        non_historical_origin_variable = True
        gradient_process_non_hist_non_hist = KratosMultiphysics.ComputeNonHistoricalNodalGradientProcess(main_model_part,
        KratosMultiphysics.DISTANCE,
        KratosMultiphysics.DISTANCE_GRADIENT,
        KratosMultiphysics.NODAL_AREA,
        non_historical_origin_variable)
        gradient_process_non_hist_non_hist.Execute()

        for node in main_model_part.Nodes:
            distance_gradient = node.GetValue(KratosMultiphysics.DISTANCE_GRADIENT)
            ref_gradient = reference_values[str(node.Id)]
            for gradient_i, gradient_i_ref in zip(distance_gradient, ref_gradient):
                self.assertAlmostEqual(gradient_i, gradient_i_ref)

        gradient_parameters = KratosMultiphysics.Parameters(""" {
            "origin_variable"                : "DISTANCE",
            "gradient_variable"              : "DISTANCE_GRADIENT",
            "non_historical_origin_variable" :  false
        }
        """)
        gradient_process = KratosMultiphysics.ComputeNonHistoricalNodalGradientProcess(main_model_part,
        gradient_parameters)
        gradient_process.Execute()

        for node in main_model_part.Nodes:
            distance_gradient = node.GetValue(KratosMultiphysics.DISTANCE_GRADIENT)
            ref_gradient = reference_values[str(node.Id)]
            for i_gradient, i_reference in zip(distance_gradient, ref_gradient):
                self.assertAlmostEqual(i_gradient, i_reference)

        gradient_parameters = KratosMultiphysics.Parameters(""" {
            "origin_variable"                : "DISTANCE",
            "gradient_variable"              : "DISTANCE_GRADIENT",
            "non_historical_origin_variable" :  true
        }
        """)
        gradient_process = KratosMultiphysics.ComputeNonHistoricalNodalGradientProcess(main_model_part,
        gradient_parameters)
        gradient_process.Execute()

        for node in main_model_part.Nodes:
            distance_gradient = node.GetValue(KratosMultiphysics.DISTANCE_GRADIENT)
            ref_gradient = reference_values[str(node.Id)]
            for i_gradient, i_reference in zip(distance_gradient, ref_gradient):
                self.assertAlmostEqual(i_gradient, i_reference)

if __name__ == '__main__':
    KratosUnittest.main()
