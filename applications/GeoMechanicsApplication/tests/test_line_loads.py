import os

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.StructuralMechanicsApplication as KratosStructural
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

    def test_nonuniform_line_load(self):
        test_name = 'non-uniform_line_load'
        parent_name = 'line_load_tests'
        file_path = test_helper.get_file_path(os.path.join(parent_name, test_name))
        simulation = test_helper.run_kratos(file_path)

        bottom_node_ids = [1, 2, 6, 11, 17, 25, 34, 46, 59, 75, 90]
        reactions_y = [reaction[1] for reaction in test_helper.get_nodal_variable(simulation, Kratos.REACTION, bottom_node_ids)]

        self.assertAlmostEqual(20.0, sum(reactions_y), 4)

        line_load_y_by_node_id = [( 91,   0.0),
                                  ( 93,   0.0),
                                  ( 98,  -4.0),
                                  (102,  -9.0),
                                  (106, -10.0),
                                  (111, -10.0),
                                  (114,  -6.0),
                                  (118,  -1.0),
                                  (122,   0.0),
                                  (124,   0.0),
                                  (125,   0.0)]
        top_node_ids = [item[0] for item in line_load_y_by_node_id]
        line_loads_y = [item[1] for item in test_helper.get_nodal_variable(simulation, KratosStructural.LINE_LOAD, top_node_ids)]
        for id_and_load_value, output_load_value in zip(line_load_y_by_node_id, line_loads_y):
            self.assertAlmostEqual(id_and_load_value[1], output_load_value, 4)

if __name__ == '__main__':
    KratosUnittest.main()
