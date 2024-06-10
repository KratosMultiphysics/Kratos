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
                                                                             node_ids=[1, 5])
        expected_heads    = [ 1., 0. ]
        for head, expected_head in zip(heads, expected_heads):
            self.assertAlmostEqual(head, expected_head)
        
if __name__=="__main__":
    KratosUnittest.main()
