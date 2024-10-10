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
        input_filename = GetFilePath("test_files/mdpa_files/test_mpi_serializer")
        ReadModelPart(input_filename, main_model_part)

        for node in main_model_part.Nodes:
            distance = node.X**2+node.Y**2+node.Z**2 - 1
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE,distance)
            node.SetValue(KratosMultiphysics.DISTANCE,distance)

        ## Read reference
        file_name = "test_files/reference_files/test_compute_nodal_gradient_process_results.json"
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
        input_filename = GetFilePath("test_files/mdpa_files/test_mpi_serializer")
        ReadModelPart(input_filename, main_model_part)

        for node in main_model_part.Nodes:
            distance = node.X**2+node.Y**2+node.Z**2 - 1
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE,distance)
            node.SetValue(KratosMultiphysics.DISTANCE,distance)

        ## Read reference
        file_name = "test_files/reference_files/test_compute_nodal_gradient_process_results.json"
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

    def testDiscontinuousDistanceProcessCutOnEdge2D(self):

        communicator = KratosMultiphysics.Testing.GetDefaultDataCommunicator()
        model = KratosMultiphysics.Model()
        main_model_part = model.CreateModelPart('main')
        main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]=2
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)

        input_filename = GetFilePath("test_files/mdpa_files/test_mpi_distance")
        ReadModelPart(input_filename, main_model_part)
        skin_model_part  = model.CreateModelPart('skin')
        skin_model_part.CreateNewNode(1, -0.25, 0.0, 0.0)
        skin_model_part.CreateNewNode(2,  0.4, 0.0, 0.0)
        skin_model_part.CreateNewElement("Element2D2N", 1, [1,2], KratosMultiphysics.Properties(0))

        params = KratosMultiphysics.Parameters()
        params.AddBool("calculate_elemental_edge_distances", True)

        KratosMultiphysics.CalculateDiscontinuousDistanceToSkinProcess2D(
            main_model_part,
            skin_model_part,
            params).Execute()

        n_intersected, n_incised = _GetNumberIncisedAndIntersectedElements(main_model_part, 3)

        n_intersected = communicator.SumAll(n_intersected)
        n_incised = communicator.SumAll(n_incised)
        self.assertEqual(n_intersected, 1)
        self.assertEqual(n_incised, 2)

        for elem in main_model_part.Elements:
            edge_distances = elem.GetValue(KratosMultiphysics.ELEMENTAL_EDGE_DISTANCES)
            if elem.Id == 2:
                self.assertVectorAlmostEqual(edge_distances, [-1.0, -1.0, -1.0])
            if elem.Id == 5:
                self.assertVectorAlmostEqual(edge_distances, [-1.0, -1.0, -1.0])

    def testDiscontinuousDistanceProcessCutThroughNode2D(self):

        communicator = KratosMultiphysics.Testing.GetDefaultDataCommunicator()
        model = KratosMultiphysics.Model()
        skin_model_part  = model.CreateModelPart('skin')
        skin_model_part.CreateNewNode(1, -0.4,  0.2, 0.0)
        skin_model_part.CreateNewNode(2,  0.4, -0.2, 0.0)
        skin_model_part.CreateNewElement("Element2D2N", 1, [1,2], KratosMultiphysics.Properties(0));

        main_model_part = model.CreateModelPart('main')
        main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]=2
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)

        input_filename = GetFilePath("test_files/mdpa_files/test_mpi_distance")
        ReadModelPart(input_filename, main_model_part)

        params = KratosMultiphysics.Parameters()
        params.AddBool("calculate_elemental_edge_distances", True)
        params.AddBool("calculate_elemental_edge_distances_extrapolated", True)
        KratosMultiphysics.CalculateDiscontinuousDistanceToSkinProcess2D(
            main_model_part,
            skin_model_part,
            params).Execute()

        n_intersected, n_incised = _GetNumberIncisedAndIntersectedElements(main_model_part, 3)

        n_intersected = communicator.SumAll(n_intersected)
        n_incised = communicator.SumAll(n_incised)
        self.assertEqual(n_intersected, 4)
        self.assertEqual(n_incised, 2)

        for elem in main_model_part.Elements:
            elem_distances = elem.GetValue(KratosMultiphysics.ELEMENTAL_DISTANCES)
            edge_distances = elem.GetValue(KratosMultiphysics.ELEMENTAL_EDGE_DISTANCES)
            edge_distances_extrapolated = elem.GetValue(KratosMultiphysics.ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED)
            if elem.Id == 2:
                self.assertVectorAlmostEqual(edge_distances, [-1.0,-1.0, -1.0])
            if elem.Id == 3:
                self.assertVectorAlmostEqual(edge_distances, [-1.0, -1.0, 1.0/3.0])
                self.assertVectorAlmostEqual(edge_distances_extrapolated, [-1.0,0.5,-1.0])
            if elem.Id == 4:
                self.assertVectorAlmostEqual(elem_distances, [0.22360679,0.0,-0.447213595])
                self.assertVectorAlmostEqual(edge_distances, [0.0, 2.0 / 3.0, -1.0])
                self.assertVectorAlmostEqual(edge_distances_extrapolated, [-1.0, -1.0, -1.0])
            if elem.Id == 6:
                self.assertVectorAlmostEqual(elem_distances, [0.447213595,0.22360679,-0.22360679])
                self.assertVectorAlmostEqual(edge_distances, [-1.0, 1.0 / 3.0 , -1.0])
                self.assertVectorAlmostEqual(edge_distances_extrapolated, [0.5, -1.0, -1.0])

def _GetNumberIncisedAndIntersectedElements(main_model_part, n_edges):
    n_intersected = 0
    n_incised = 0
    for elem in main_model_part.Elements:
        n_cut_edges = 0
        edge_dist = elem.GetValue(KratosMultiphysics.ELEMENTAL_EDGE_DISTANCES)
        for distance in edge_dist:
            if distance >= 0.0:
                n_cut_edges += 1
        if n_cut_edges>0:
            if n_cut_edges>1:
                n_intersected += 1
            else:
                n_incised += 1
    return n_intersected, n_incised

if __name__ == '__main__':
    KratosUnittest.main()
