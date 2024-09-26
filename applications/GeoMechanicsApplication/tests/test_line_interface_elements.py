import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper

class KratosGeoMechanicsInterfaceElementTests(KratosUnittest.TestCase):
    """
    This class contains one test case for a 3+3 line interface element which is loaded
    in compression as well as in shear. This is done both by applying a load (Neumann)
    and prescribing a displacement (Dirichlet) on one side of the interface.
    """
    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass
        
    def test_3_plus_3_line_interface_element_with_neumann_conditions(self):
        simulation = self.run_simulation('Neumann')
        self.assert_outputs_for_3_plus_3_line_interface_element(simulation)

    def test_3_plus_3_line_interface_element_with_dirichlet_conditions(self):
        simulation = self.run_simulation('Dirichlet')
        self.assert_outputs_for_3_plus_3_line_interface_element(simulation)


    @staticmethod
    def run_simulation(condition_type):
        file_path = test_helper.get_file_path(os.path.join('line_interface_elements', f'{condition_type}_single_stage'))
        return test_helper.run_kratos(file_path)


    def assert_outputs_for_3_plus_3_line_interface_element(self, simulation):
        displacements = test_helper.get_displacement(simulation)

        # Check the horizontal element
        shear_traction = 667.0
        shear_stiffness = 1.5e7
        expected_horizontal_displacement = -shear_traction / shear_stiffness
        normal_traction = 333.0
        normal_stiffness = 3.0e7
        expected_vertical_displacement = -normal_traction / normal_stiffness
        top_node_indices = [3, 4, 5] # These correspond to nodes 11, 12, 13
        for index in top_node_indices:
            self.assertAlmostEqual(displacements[index][0], expected_horizontal_displacement)
            self.assertAlmostEqual(displacements[index][1], expected_vertical_displacement)

        tractions = test_helper.get_cauchy_stress_vectors(simulation)
        tractions_horizontal_element = tractions[0]
        for index in range(3):
            self.assertAlmostEqual(tractions_horizontal_element[index][0], -normal_traction)
            self.assertAlmostEqual(tractions_horizontal_element[index][1], -shear_traction)

        relative_displacements = test_helper.get_strain_vectors(simulation)
        relative_displacements_horizontal_element = relative_displacements[0]
        for index in range(3):
            # Notice the first index is the normal displacement (which is vertical for this horizontally oriented
            # element) and the second index is the shear displacement (which is horizontal for this element)
            self.assertAlmostEqual(relative_displacements_horizontal_element[index][0], expected_vertical_displacement)
            self.assertAlmostEqual(relative_displacements_horizontal_element[index][1], expected_horizontal_displacement)

        # Check the vertical element
        left_node_indices = [6, 7, 8] # These correspond to nodes 21, 22, 23
        expected_horizontal_displacement = normal_traction / normal_stiffness
        expected_vertical_displacement = -shear_traction / shear_stiffness
        for index in left_node_indices:
            self.assertAlmostEqual(displacements[index][0], expected_horizontal_displacement)
            self.assertAlmostEqual(displacements[index][1], expected_vertical_displacement)

        tractions_vertical_element = tractions[1]
        for index in range(3):
            self.assertAlmostEqual(tractions_vertical_element[index][0], -normal_traction)
            self.assertAlmostEqual(tractions_vertical_element[index][1], -shear_traction)

        relative_displacement_vertical_element = relative_displacements[1]
        for index in range(3):
            self.assertAlmostEqual(relative_displacement_vertical_element[index][0], -expected_horizontal_displacement)
            self.assertAlmostEqual(relative_displacement_vertical_element[index][1], expected_vertical_displacement)



if __name__ == '__main__':
    KratosUnittest.main()
