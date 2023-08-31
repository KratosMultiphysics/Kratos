import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper

class KratosGeoMechanicsLineLoadTests(KratosUnittest.TestCase):
    """
    This class contains benchmark tests to test line loads (in 3D). 
    No analytical solution is defined. This is to check if find neighbour elements work when line elements are applied in 3D.
    """

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def test_line_load_3D2N_hex(self):
        test_name = 'line_load_3D2N_hex'
        parent_name = 'line_load_tests'
        file_path = test_helper.get_file_path(os.path.join(parent_name, test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        displacements = test_helper.get_displacement(simulation)
        y_displacements = [displacement[1] for displacement in displacements]

        node_number = 2
        expected_value = 0.00042943
        self.assertAlmostEqual(expected_value, y_displacements[node_number-1], 6)

    def test_line_load_3D2N_tet(self):
        test_name = 'line_load_3D2N_tet'
        parent_name = 'line_load_tests'
        file_path = test_helper.get_file_path(os.path.join(parent_name, test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        displacements = test_helper.get_displacement(simulation)
        y_displacements = [displacement[1] for displacement in displacements]

        node_number = 4
        expected_value = 0.00034701
        self.assertAlmostEqual(expected_value, y_displacements[node_number-1], 6)

    def test_line_load_3D3N_hex(self):
        test_name = 'line_load_3D3N_hex'
        parent_name = 'line_load_tests'
        file_path = test_helper.get_file_path(os.path.join(parent_name, test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        displacements = test_helper.get_displacement(simulation)
        y_displacements = [displacement[1] for displacement in displacements]

        node_number = 5
        expected_value = 0.0006958
        self.assertAlmostEqual(expected_value, y_displacements[node_number-1], 6)

    def test_line_load_3D3N_tet(self):
        test_name = 'line_load_3D3N_tet'
        parent_name = 'line_load_tests'
        file_path = test_helper.get_file_path(os.path.join(parent_name, test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        displacements = test_helper.get_displacement(simulation)
        y_displacements = [displacement[1] for displacement in displacements]

        node_number = 11
        expected_value = 0.00064506
        self.assertAlmostEqual(expected_value, y_displacements[node_number-1], 6)

if __name__ == '__main__':
    KratosUnittest.main()