import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import GiDOutputFileReader
import test_helper


class KratosGeoMechanicsExtrapolationTests(KratosUnittest.TestCase):
    """
    Class that contains geo element integration points to nodes extrapolation tests.
    """

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def test_3_noded_triangle(self):
        directory  = 'test_integration_node_extrapolation'
        test_name  = '3_noded_triangle'
        file_path  = test_helper.get_file_path(os.path.join(directory, test_name))
        test_helper.run_kratos(file_path)
        reader     = GiDOutputFileReader()
        simulation_output = reader.read_output_from(os.path.join(file_path, test_name+'.post.res'))
        heads = GiDOutputFileReader.nodal_values_at_time(
            "HYDRAULIC_HEAD", 1, simulation_output, node_ids=[1, 3, 5]
        )
        expected_heads    = [ 1., 0.5, 0. ]
        for head, expected_head in zip(heads, expected_heads):
            self.assertAlmostEqual(head, expected_head, places=4)

    def test_4_noded_quadrilateral(self):
        directory  = 'test_integration_node_extrapolation'
        test_name = '4_noded_quadrilateral'
        file_path = test_helper.get_file_path(os.path.join(directory, test_name))
        test_helper.run_kratos(file_path)
        reader     = GiDOutputFileReader()
        simulation_output = reader.read_output_from(os.path.join(file_path, test_name+'.post.res'))
        heads = GiDOutputFileReader.nodal_values_at_time(
            "HYDRAULIC_HEAD", 1, simulation_output, node_ids=[2, 4]
        )
        expected_heads    = [ 0., 1. ]
        for head, expected_head in zip(heads, expected_heads):
            self.assertAlmostEqual(head, expected_head, places=4)

    def test_6_noded_triangle(self):
        directory  = 'test_integration_node_extrapolation'
        test_name = '6_noded_triangle'
        file_path = test_helper.get_file_path(os.path.join(directory, test_name))
        test_helper.run_kratos(file_path)
        reader     = GiDOutputFileReader()
        simulation_output = reader.read_output_from(os.path.join(file_path, test_name+'.post.res'))
        heads = GiDOutputFileReader.nodal_values_at_time(
            "HYDRAULIC_HEAD", 1, simulation_output, node_ids=[1, 3, 5]
        )
        expected_heads    = [ 1., 0.5, 0. ]
        for head, expected_head in zip(heads, expected_heads):
            self.assertAlmostEqual(head, expected_head, places=4)

    def test_8_noded_quadrilateral(self):
        directory  = 'test_integration_node_extrapolation'
        test_name = '8_noded_quadrilateral'
        file_path = test_helper.get_file_path(os.path.join(directory, test_name))
        test_helper.run_kratos(file_path)
        reader     = GiDOutputFileReader()
        simulation_output = reader.read_output_from(os.path.join(file_path, test_name+'.post.res'))
        heads = GiDOutputFileReader.nodal_values_at_time(
            "HYDRAULIC_HEAD", 1, simulation_output, node_ids=[2, 4, 8]
        )
        expected_heads    = [ 1.09636, 0.292984, 1.29298]
        for head, expected_head in zip(heads, expected_heads):
            self.assertAlmostEqual(head, expected_head, places=4)

        cauchy_stress_tensors = GiDOutputFileReader.nodal_values_at_time(
            "CAUCHY_STRESS_TENSOR", 1, simulation_output, node_ids=[2, 4, 8]
        )
        expected_tensors    = [ [-1501.6, -6006.41, -1501.6, -376.767, 0, 0],
                                [733.933, 2935.73, 733.933, -370.876, 0, 0],
                                [-8.04944, -32.1978, -8.04944, -373.822, 0, 0]]
        for cauchy_stress_tensor, expected_tensor in zip(cauchy_stress_tensors, expected_tensors):
            for cauchy_stress, expected in zip(cauchy_stress_tensor, expected_tensor):
                self.assertAlmostEqual(cauchy_stress, expected, places=None, delta=0.01)

        fluid_flux_vectors = GiDOutputFileReader.nodal_values_at_time(
            "FLUID_FLUX_VECTOR", 1, simulation_output, node_ids=[2, 4, 8]
        )
        expected_fluid_flux_vectors    = [ [1.38449e-05, -4.52218e-07, 0],
                            [-3.51008e-06, 2.17395e-05, 0],
                            [5.16742e-06, 1.157e-05, 0]]
        for fluid_flux_vector, expected_fluid_flux_vector in zip(fluid_flux_vectors, expected_fluid_flux_vectors):
            for fluid_flux_component, expected in zip(fluid_flux_vector, expected_fluid_flux_vector):
                self.assertAlmostEqual(fluid_flux_component, expected, places=4)

if __name__=="__main__":
    KratosUnittest.main()
