import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import GiDOutputFileReader
import KratosMultiphysics.GeoMechanicsApplication.run_multiple_stages as run_multiple_stages
import test_helper

class KratosGeoMechanicsLabElementTests(KratosUnittest.TestCase):
    """
    This class contains some element tests, such as triaxial and oedometer tests
    """
    def test_triaxial(self):
        """
        Regression test for the triaxial experiment.
        """
        expected_disp = [[0.0, -0.2, 0.0], [0.0527776, -0.2, 0.0], [0.0, -0.100033, 0.0], [0.0524025, -0.0996931, 0.0], [0.0, 0.0, 0.0], [0.105197, -0.2, 0.0], [0.105114, -0.100049, 0.0], [0.0524406, 0.0, 0.0], [0.104632, 0.0, 0.0]]
        expected_stress = [[-99.9808, -252.622, -99.9806, 0.193199, 0.0, 0.0], [-99.9991, -252.668, -99.9991, 0.00846584, 0, 0]]
        expected_strain = [[0.104863, -0.19973, 0.104946, 0.000440186, 0.0, 0.0], [0.1055, -0.200303, 0.104922, 3.84218e-05, 0, 0]]
        self._run_triaxial_regression_test('drained', 'triaxial_test_output.post.res', expected_disp, expected_stress,
                                           expected_strain, 4)


    def test_triaxial_undrained(self):
        """
        Regression test for the undrained triaxial experiment.
        """
        expected_disp = [[0.0, -0.2, 0.0], [0.045, -0.2, 0.0], [0.0, -0.1, 0.0], [ 0.045, -0.1, 0.0],[0.0, 0.0, 0.0], [0.09, -0.2, 0.0], [0.09, -0.1, 0.0], [0.045, 0.0, 0.0], [0.09, 0.0, 0.0]]
        expected_stress = [[-100.0, -4740.0, -100.0, 0.0, 0.0, 0.0], [-100.0, -4740.0, -100.0, 0.0, 0.0, 0.0]]
        expected_strain = [[0.09, -0.2, 0.09, 0.0, 0.0, 0.0], [0.09, -0.2, 0.09, 0.0, 0.0, 0.0]]
        self._run_triaxial_regression_test('undrained', 'triaxial_undrained_test_output.post.res', expected_disp,
                                           expected_stress, expected_strain, 4)

    def _run_triaxial_regression_test(self, stage_name, output_file_name, expected_disp, expected_stress,
                                      expected_strain, precision_places):
        test_name = os.path.join('test_triaxial', stage_name)
        file_path = test_helper.get_file_path(os.path.join('test_element_lab', test_name))
        test_helper.run_kratos(file_path)

        # Read the output files from the simulation for comparison.
        reader = GiDOutputFileReader()
        result = reader.read_output_from(os.path.join(file_path, output_file_name))
        time = 1.0

        # Assert the displacement in all nodes in all directions.
        for node_id, expected_node_displacement in enumerate(expected_disp, start=1):
            node_displacement = reader.nodal_values_at_time("DISPLACEMENT", time, result, [node_id])[0]
            self.assertVectorAlmostEqual(node_displacement, expected_node_displacement, precision_places)

        self._assert_first_ip_tensor(reader, result, "CAUCHY_STRESS_TENSOR", expected_stress, precision_places, time)
        self._assert_first_ip_tensor(reader, result, "ENGINEERING_STRAIN_TENSOR", expected_strain, precision_places, time)

    def _assert_first_ip_tensor(self, reader, result, variable_name, expected_values, precision_places, time=1.0):
        for element_id, expected_tensor in enumerate(expected_values, start=1):
            tensor = reader.element_integration_point_values_at_time(variable_name, time, result, [element_id], [0])[0]
            self._assert_integration_point_tensor_results(tensor, expected_tensor, precision_places, variable_name)

    def test_triaxial_comp_6n(self):
        """
        Drained compression triaxial test on Mohr-Coulomb model with axisymmetric 2D6N elements
        It consistes of two calculation phases:
        1) apply confining stress of -100 kPa
        2) apply deviatoric stress of -200 kPa
        """
        project_path = test_helper.get_file_path(os.path.join('test_element_lab', 'triaxial_comp_6n'))

        n_stages = 2
        run_multiple_stages.run_stages(project_path, n_stages)

        reader = GiDOutputFileReader()

        # Assert
        output_data = reader.read_output_from(os.path.join(project_path, "triaxial_comp_6n_stage1.post.res"))
        time = 1.0
        stress_vectors_per_element = reader.element_integration_point_values_at_time("CAUCHY_STRESS_TENSOR", time, output_data)
        number_of_elements = 2
        self.assertEqual(number_of_elements, len(stress_vectors_per_element))

        number_of_integration_points_per_element = 3
        expected_stress_vector = [-100.0, -100.0, -100.0]
        for element_stress_vectors in stress_vectors_per_element:
            self.assertEqual(number_of_integration_points_per_element, len(element_stress_vectors))
            self._assert_integration_point_tensor_results(element_stress_vectors, expected_stress_vector, 3,"CAUCHY_STRESS_TENSOR")

        output_data = reader.read_output_from(os.path.join(project_path, "triaxial_comp_6n_stage2.post.res"))
        time = 1.25
        stress_vectors_per_element = reader.element_integration_point_values_at_time("CAUCHY_STRESS_TENSOR", time, output_data)
        self.assertEqual(number_of_elements, len(stress_vectors_per_element))

        expected_stress_vector = [-100.0, -300.0, -100.0]
        for element_stress_vectors in stress_vectors_per_element:
            self.assertEqual(number_of_integration_points_per_element, len(element_stress_vectors))
            self._assert_integration_point_tensor_results(element_stress_vectors, expected_stress_vector, 2,"CAUCHY_STRESS_TENSOR")


    def test_oedometer_ULFEM(self):
        """
        Oedometer test on a linear elastic model with 2D6N elements
        """
        test_name = 'oedometer_ULFEM'
        project_path = test_helper.get_file_path(os.path.join('test_element_lab', test_name))
        simulation = test_helper.run_kratos(project_path)
        effective_stresses = test_helper.get_cauchy_stress_tensor(simulation)
        effective_stresses_xx = [integration_point[0,0] for element in effective_stresses for integration_point in element]
        effective_stresses_yy = [integration_point[1,1] for element in effective_stresses for integration_point in element]
        effective_stresses_zz = [integration_point[2,2] for element in effective_stresses for integration_point in element]

        # Assert integration point information
        for idx, effective_stress_xx in enumerate(effective_stresses_xx):
            self.assertAlmostEqual(0.0,      effective_stress_xx,        3)
            self.assertAlmostEqual(-1000000, effective_stresses_yy[idx], 2)
            self.assertAlmostEqual(0.0,      effective_stresses_zz[idx], 3)

        top_node_nbrs = [1, 2]
        output_file_path = os.path.join(project_path, test_name+'.post.res')
        output_reader = GiDOutputFileReader()
        output_data = output_reader.read_output_from(output_file_path)
        displacements = GiDOutputFileReader.nodal_values_at_time("DISPLACEMENT", 0.1, output_data, node_ids=top_node_nbrs)
        for displacement in displacements:
            y_displacement = displacement[1]
            self.assertAlmostEqual(-0.00990099, y_displacement, 6)

        displacements = GiDOutputFileReader.nodal_values_at_time("DISPLACEMENT", 0.7, output_data, node_ids=top_node_nbrs)
        for displacement in displacements:
            y_displacement = displacement[1]
            self.assertAlmostEqual(-0.0654206, y_displacement, 6)

        displacements = GiDOutputFileReader.nodal_values_at_time("DISPLACEMENT", 1.0, output_data, node_ids=top_node_nbrs)
        for displacement in displacements:
            y_displacement = displacement[1]
            self.assertAlmostEqual(-0.0909090909516868, y_displacement, 6)

    def test_oedometer_ULFEM_diff_order(self):
        """
        Oedometer test on a linear elastic model with 2D6N with different order elements
        """
        test_name = 'oedometer_ULFEM_diff_order'
        project_path = test_helper.get_file_path(os.path.join('test_element_lab', test_name))
        simulation = test_helper.run_kratos(project_path)
        effective_stresses = test_helper.get_cauchy_stress_tensor(simulation)
        displacements = test_helper.get_displacement(simulation)

        effective_stresses_xx = [integration_point[0,0] for element in effective_stresses for integration_point in element]
        effective_stresses_yy = [integration_point[1,1] for element in effective_stresses for integration_point in element]
        effective_stresses_zz = [integration_point[2,2] for element in effective_stresses for integration_point in element]

        # Assert integration point information
        for idx, effective_stress_xx in enumerate(effective_stresses_xx):
            self.assertAlmostEqual(0.0,   effective_stress_xx,        3)
            self.assertAlmostEqual(-1000, effective_stresses_yy[idx], 3)
            self.assertAlmostEqual(0.0,   effective_stresses_zz[idx], 3)

        y_displacements = [displacement[1] for displacement in displacements]
        top_node_nbrs = [1]
        for top_node_nbr in top_node_nbrs:
            self.assertAlmostEqual(-1e-04, y_displacements[top_node_nbr], 6)

    def _assert_integration_point_tensor_results(self, integration_point_tensors, expected_integration_point_tensor, places, result_name):
        for idx, ip_tensor in enumerate(integration_point_tensors):
            self.assertVectorAlmostEqual(expected_integration_point_tensor[:2], ip_tensor[:2], places, msg = f"{result_name} component xx, yy, zz at integration point {idx}")


if __name__ == '__main__':
    KratosUnittest.main()
