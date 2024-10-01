import os

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper

class KratosGeoMechanicsInterfaceElementTests(KratosUnittest.TestCase):
    """
    This class contains one test case for a 3+3 line interface element which is loaded
    in compression as well as in shear. This is done both by applying a load (Neumann)
    and prescribing a displacement (Dirichlet) on one side of the interface.
    """
    def setUp(self):
        self.model = Kratos.Model()
        self.normal_stiffness = 3.0e7
        self.shear_stiffness = 1.5e7

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
        expected_horizontal_displacement = -shear_traction / self.shear_stiffness
        normal_traction = 333.0
        expected_vertical_displacement = -normal_traction / self.normal_stiffness
        top_node_indices = [3, 4, 5] # These correspond to nodes 11, 12, 13
        for index in top_node_indices:
            self.assertAlmostEqual(displacements[index][0], expected_horizontal_displacement)
            self.assertAlmostEqual(displacements[index][1], expected_vertical_displacement)

        tractions = test_helper.get_on_integration_points(simulation, Kratos.CAUCHY_STRESS_VECTOR)
        tractions_horizontal_element = tractions[0]
        for index in range(3):
            self.assertAlmostEqual(tractions_horizontal_element[index][0], -normal_traction)
            self.assertAlmostEqual(tractions_horizontal_element[index][1], -shear_traction)

        relative_displacements = test_helper.get_on_integration_points(simulation, Kratos.STRAIN)
        relative_displacements_horizontal_element = relative_displacements[0]
        for index in range(3):
            # Notice the first index is the normal displacement (which is vertical for this horizontally oriented
            # element) and the second index is the shear displacement (which is horizontal for this element)
            self.assertAlmostEqual(relative_displacements_horizontal_element[index][0], expected_vertical_displacement)
            self.assertAlmostEqual(relative_displacements_horizontal_element[index][1], expected_horizontal_displacement)

        # Check the vertical element
        left_node_indices = [6, 7, 8] # These correspond to nodes 21, 22, 23
        expected_horizontal_displacement = normal_traction / self.normal_stiffness
        expected_vertical_displacement = -shear_traction / self.shear_stiffness
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


    def assertVectorsAlmostEqual(self, actual_vectors, expected_vectors):
        for actual_vector, expected_vector in zip(actual_vectors, expected_vectors):
            self.assertVectorAlmostEqual(actual_vector, expected_vector)


    def assert_results_of_multi_stage_test(self, stage, displacement_vector):
        expected_displacement_vectors = [[0.0, 0.0, 0.0]] * 3  # the first three nodes have been fixed
        expected_displacement_vectors += [displacement_vector] * 3  # the last three nodes have prescribed non-zero displacements
        self.assertVectorsAlmostEqual(test_helper.get_displacement(stage), expected_displacement_vectors)

        expected_normal_relative_displacement = displacement_vector[1]
        expected_tangential_relative_displacement = displacement_vector[0]
        expected_relative_displacement_vectors = [[expected_normal_relative_displacement, expected_tangential_relative_displacement]] * 3
        self.assertVectorsAlmostEqual(test_helper.get_on_integration_points(stage, Kratos.STRAIN)[0],
                                      expected_relative_displacement_vectors)

        expected_traction_vectors = [[self.normal_stiffness * expected_normal_relative_displacement,
                                      self.shear_stiffness * expected_tangential_relative_displacement]] * 3
        self.assertVectorsAlmostEqual(test_helper.get_on_integration_points(stage, Kratos.CAUCHY_STRESS_VECTOR)[0],
                                      expected_traction_vectors)


    def test_multi_stage_3_plus_3_line_interface_element_with_dirichlet_conditions(self):
        file_path = test_helper.get_file_path(os.path.join('line_interface_elements', 'Dirichlet_multi_stage'))

        initial_cwd = os.getcwd()
        os.chdir(file_path)

        project_parameters_file_names = ['ProjectParameters_stage1.json', 'ProjectParameters_stage2.json']
        prescribed_displacement_vectors = [[-8.8933333333333332e-5, -2.22e-5, 0.0], [-2.0e-4, -4.44e-5, 0.0]]
        for file_name, displacement_vector in zip(project_parameters_file_names, prescribed_displacement_vectors):
            stage = test_helper.make_geomechanics_analysis(self.model, os.path.join(file_path, file_name))
            stage.Run()

            self.assert_results_of_multi_stage_test(stage, displacement_vector)

        os.chdir(initial_cwd)


    def test_multi_stage_3_plus_3_line_interface_element_with_neumann_conditions(self):
        file_path = test_helper.get_file_path(os.path.join('line_interface_elements', 'Neumann_multi_stage'))

        initial_cwd = os.getcwd()
        os.chdir(file_path)

        project_parameters_file_names = ['ProjectParameters_stage1.json', 'ProjectParameters_stage2.json']
        expected_displacement_vectors_of_loaded_side = [[-8.8933333333333332e-5, -2.22e-5, 0.0], [-2.0e-4, -4.44e-5, 0.0]]
        for file_name, displacement_vector in zip(project_parameters_file_names, expected_displacement_vectors_of_loaded_side):
            stage = test_helper.make_geomechanics_analysis(self.model, os.path.join(file_path, file_name))
            stage.Run()

            self.assert_results_of_multi_stage_test(stage, displacement_vector)

        os.chdir(initial_cwd)

if __name__ == '__main__':
    KratosUnittest.main()
