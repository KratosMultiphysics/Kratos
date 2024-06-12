import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
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
        reader     = test_helper.GiDOutputFileReader()
        simulation_output = reader.read_output_from(os.path.join(file_path, test_name+'.post.res'))
        heads             = test_helper.GiDOutputFileReader.nodal_values_at_time("HYDRAULIC_HEAD", 1, simulation_output,
                                                                             node_ids=[1, 3, 5])
        expected_heads    = [ 1., 0.5, 0. ]
        for head, expected_head in zip(heads, expected_heads):
            self.assertAlmostEqual(head, expected_head)
       
    def test_4_noded_quadrilateral(self):
        directory  = 'test_integration_node_extrapolation'
        test_name = '4_noded_quadrilateral'
        file_path = test_helper.get_file_path(os.path.join(directory, test_name))
        test_helper.run_kratos(file_path)
        reader     = test_helper.GiDOutputFileReader()
        simulation_output = reader.read_output_from(os.path.join(file_path, test_name+'.post.res'))
        heads             = test_helper.GiDOutputFileReader.nodal_values_at_time("HYDRAULIC_HEAD", 1, simulation_output,
                                                                             node_ids=[1, 5])
        expected_heads    = [ 1., 0. ]
        for head, expected_head in zip(heads, expected_heads):
            self.assertAlmostEqual(head, expected_head)

    def test_6_noded_triangle(self):
        directory  = 'test_integration_node_extrapolation'
        test_name = '6_noded_triangle'
        file_path = test_helper.get_file_path(os.path.join(directory, test_name))
        test_helper.run_kratos(file_path)
        reader     = test_helper.GiDOutputFileReader()
        simulation_output = reader.read_output_from(os.path.join(file_path, test_name+'.post.res'))
        heads             = test_helper.GiDOutputFileReader.nodal_values_at_time("HYDRAULIC_HEAD", 1, simulation_output,
                                                                             node_ids=[1, 3, 5])
        expected_heads    = [ 1., 0.5, 0. ]
        for head, expected_head in zip(heads, expected_heads):
            self.assertAlmostEqual(head, expected_head)

    def test_8_noded_quadrilateral(self):
        directory  = 'test_integration_node_extrapolation'
        test_name = '8_noded_quadrilateral'
        file_path = test_helper.get_file_path(os.path.join(directory, test_name))
        test_helper.run_kratos(file_path)
        reader     = test_helper.GiDOutputFileReader()
        simulation_output = reader.read_output_from(os.path.join(file_path, test_name+'.post.res'))
        heads             = test_helper.GiDOutputFileReader.nodal_values_at_time("HYDRAULIC_HEAD", 1, simulation_output,
                                                                             node_ids=[1, 5, 9])
        expected_heads    = [ -0.396837, -1.39684, -0.896837]
        for head, expected_head in zip(heads, expected_heads):
            self.assertAlmostEqual(head, expected_head)


        cauchy_stress_tensors             = test_helper.GiDOutputFileReader.nodal_values_at_time("CAUCHY_STRESS_TENSOR", 1, simulation_output,
                                                                                     node_ids=[1, 5, 9])
        expected_tensors    = [ [10915.7, 2728.92, 2728.92, 4094.93, 0, 0],
                                [-10899.1, -2724.78, -2724.78, 8185.21, 0, 0],
                                [8.28365, 2.07091, 2.07091, 2045.14, 0, 0]]
        for cauchy_stress_tensor, expected_tensor in zip(cauchy_stress_tensors, expected_tensors):
            for (cauchy_stress, expected) in zip(cauchy_stress_tensor, expected_tensor):
                self.assertAlmostEqual(cauchy_stress, expected)


        fluid_flux_vectors             = test_helper.GiDOutputFileReader.nodal_values_at_time("FLUID_FLUX_VECTOR", 1, simulation_output,
                                                                                                         node_ids=[1, 5, 9])
        expected_fluid_flux_vectors    = [ [-5.26143e-22, 5.42692e-05, 0],
                            [-5.26143e-22, -4.26992e-05, 0],
                            [-5.26143e-22, 5.785e-06, 0]]
        for fluid_flux_vector, expected_fluid_flux_vector in zip(fluid_flux_vectors, expected_fluid_flux_vectors):
            for (fluid_flux_component, expected) in zip(fluid_flux_vector, expected_fluid_flux_vector):
                self.assertAlmostEqual(fluid_flux_component, expected)

if __name__=="__main__":
    KratosUnittest.main()
