from csv import reader
import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import GiDOutputFileReader
import KratosMultiphysics.GeoMechanicsApplication.run_multiple_stages as run_multiple_stages
import test_helper

class KratosGeoMechanicsLabElementTests(KratosUnittest.TestCase):
    """
    This class contains some element tests, such as triaxial and oedometer tests.
    """
    def test_triaxial_drained(self):
        """Regression test for the triaxial experiment on a mohr coulomb model with a constant pore water pressure."""
        expected_disp = [[0.0, -0.2, 0.0], [0.0527776, -0.2, 0.0], [0.0, -0.100033, 0.0], [0.0524025, -0.0996931, 0.0], [0.0, 0.0, 0.0], [0.105197, -0.2, 0.0], [0.105114, -0.100049, 0.0], [0.0524406, 0.0, 0.0], [0.104632, 0.0, 0.0]]
        expected_stress = [[-99.9808, -252.622, -99.9806, 0.193199, 0.0, 0.0], [-99.8828, -252.381, -99.8827, 0.0448209, 0, 0], [-100.045, -252.78, -100.045, 0.0314238, 0, 0], [-99.9991, -252.668, -99.9991, 0.00846584, 0, 0], [-99.967, -252.589, -99.967, -0.0103576, 0, 0], [-100.064, -252.827, -100.064, 0.137162, 0, 0]]
        expected_strain = [[0.104863, -0.19973, 0.104946, 0.000440186, 0.0, 0.0], [0.104939, -0.199861, 0.105022, 0.000100528, 0, 0], [0.10444, -0.199182, 0.104773, 6.24562e-05, 0, 0], [0.1055, -0.200303, 0.104922, 3.84218e-05, 0, 0], [0.104915, -0.200106, 0.105298, -4.97694e-06, 0, 0], [0.105631, -0.200818, 0.105335, 0.000287088, 0, 0]]
        self._run_triaxial_regression_test('drained', 'triaxial_test_output.post.res', expected_disp, expected_stress,
                                           expected_strain, 4, assert_all_integration_points=True)

    def test_triaxial_undrained(self):
        """
        Regression test for the undrained triaxial experiment.
        """
        expected_disp = [[0.0, -0.2, 0.0], [0.045, -0.2, 0.0], [0.0, -0.1, 0.0], [ 0.045, -0.1, 0.0],[0.0, 0.0, 0.0], [0.09, -0.2, 0.0], [0.09, -0.1, 0.0], [0.045, 0.0, 0.0], [0.09, 0.0, 0.0]]
        expected_stress = [[-100.0, -4740.0, -100.0, 0.0, 0.0, 0.0], [-100.0, -4740.0, -100.0, 0.0, 0.0, 0.0]]
        expected_strain = [[0.09, -0.2, 0.09, 0.0, 0.0, 0.0], [0.09, -0.2, 0.09, 0.0, 0.0, 0.0]]
        self._run_triaxial_regression_test('undrained', 'triaxial_undrained_test_output.post.res', expected_disp,
                                           expected_stress, expected_strain, 4)

    def _run_triaxial_regression_test(self, stage_name, output_file_name, expected_displacement, expected_stress,
                                      expected_strain, precision_places, assert_all_integration_points=False):
        test_name = 'test_triaxial'
        file_path = test_helper.get_file_path(os.path.join('test_element_lab', test_name, stage_name))
        test_helper.run_kratos(file_path)

        # Read the output files from the simulation for comparison.
        reader = GiDOutputFileReader()
        result = reader.read_output_from(os.path.join(file_path, output_file_name))
        time = 1.0

        # Assert the displacement in all nodes in all directions.
        for node_id, expected_node_displacement in enumerate(expected_displacement, start=1):
            node_displacement = reader.nodal_values_at_time("DISPLACEMENT", time, result, [node_id])[0]
            self.assertVectorAlmostEqual(node_displacement, expected_node_displacement, precision_places)

        if assert_all_integration_points:
            self._assert_integration_point_tensors(reader, result, "CAUCHY_STRESS_TENSOR", expected_stress,
                                                  precision_places, time)
            self._assert_integration_point_tensors(reader, result, "ENGINEERING_STRAIN_TENSOR", expected_strain,
                                                  precision_places, time)
        else:
            self._assert_first_ip_tensor(reader, result, "CAUCHY_STRESS_TENSOR", expected_stress,
                                         precision_places, time)
            self._assert_first_ip_tensor(reader, result, "ENGINEERING_STRAIN_TENSOR", expected_strain,
                                         precision_places, time)
            
    def test_oedometer_drained(self):
        """Regression test for the oedometer experiment on a linear elastic model with constant pore water pressure."""
        expected_stress_per_ip = [-333.333, -1000.0, -333.333, 0.0, 0.0, 0.0]
        expected_stress = [expected_stress_per_ip] * 6
        expected_y_displacements = [-0.0208333, -0.0416667, -0.0625, -0.0833333]
        displacement_times = [0.25, 0.5, 0.75, 1.0]
        top_nodes = [7, 8, 9]
        self._run_oedometer_regression_test('drained', 'test_oedometer_output.post.res', expected_stress,
                                            expected_y_displacements, displacement_times, top_nodes, 6)
        
    def _run_oedometer_regression_test(self, stage_name, output_file_name, expected_stress,
                                       expected_y_displacements, displacement_times, top_nodes, precision_places):
        """Run the oedometer regression test and validate stress tensors and y-displacements in time."""
        test_name = 'test_oedometer'
        file_path = test_helper.get_file_path(os.path.join('test_element_lab', test_name, stage_name))
        test_helper.run_kratos(file_path)

        reader = GiDOutputFileReader()
        result = reader.read_output_from(os.path.join(file_path, output_file_name))

        self._assert_integration_point_tensors(reader, result, "CAUCHY_STRESS_TENSOR", expected_stress,
                                              precision_places, time=1.0)

        for time, expected_y in zip(displacement_times, expected_y_displacements):
            self._assert_y_displacements_at_time(reader, result, top_nodes, expected_y, precision_places, time)

    def test_triaxial_comp_6n(self):
        """
        Drained compression triaxial test on Mohr-Coulomb model with axisymmetric 2D6N elements
        It consists of two calculation phases:
        1) apply confining stress of -100 kPa
        2) apply deviatoric stress of -200 kPa
        """
        test_name = 'triaxial_comp_6n'
        project_path = test_helper.get_file_path(os.path.join('test_element_lab', test_name))
        n_stages = 2
        run_multiple_stages.run_stages(project_path, n_stages)
        reader = GiDOutputFileReader()

        result = reader.read_output_from(os.path.join(project_path, "triaxial_comp_6n_stage1.post.res"))
        expected_stress = [[-100.0, -100.0, -100.0, 0.0, 0.0, 0.0]] * 6
        self._assert_integration_point_tensors(reader, result, "CAUCHY_STRESS_TENSOR", expected_stress, 3, time=1.0)

        result = reader.read_output_from(os.path.join(project_path, "triaxial_comp_6n_stage2.post.res"))
        expected_stress = [[-100.0, -300.0, -100.0, 0.0, 0.0, 0.0]] * 6
        self._assert_integration_point_tensors(reader, result, "CAUCHY_STRESS_TENSOR", expected_stress, 2, time=1.25)

    def test_oedometer_ULFEM(self):
        """
        Oedometer test on a linear elastic model with 2D6N elements
        """
        test_name = 'oedometer_ULFEM'
        file_path = test_helper.get_file_path(os.path.join('test_element_lab', test_name))
        simulation = test_helper.run_kratos(file_path)
        effective_stresses = test_helper.get_cauchy_stress_tensor(simulation)
        self._assert_oedometer_effective_stresses(effective_stresses, -1000000.0, 2)

        top_node_nbrs = [1, 2]
        output_file = os.path.join(file_path, test_name+'.post.res')
        reader = GiDOutputFileReader()
        result = reader.read_output_from(output_file)

        self._assert_y_displacements_at_time(reader, result, top_node_nbrs, -0.00990099, 6, 0.1)
        self._assert_y_displacements_at_time(reader, result, top_node_nbrs, -0.0654206, 6, 0.7)
        self._assert_y_displacements_at_time(reader, result, top_node_nbrs, -0.0909090909516868, 6, 1.0)

    def test_oedometer_ULFEM_diff_order(self):
        """
        Oedometer test on a linear elastic model with 2D6N with different order elements
        """
        test_name = 'oedometer_ULFEM_diff_order'
        project_path = test_helper.get_file_path(os.path.join('test_element_lab', test_name))
        simulation = test_helper.run_kratos(project_path)
        effective_stresses = test_helper.get_cauchy_stress_tensor(simulation)
        self._assert_oedometer_effective_stresses(effective_stresses, -1000.0, 3)

        output_file = os.path.join(project_path, test_name + '.post.res')
        reader = GiDOutputFileReader()
        result = reader.read_output_from(output_file)
        top_node_nbrs = [1]
        self._assert_y_displacements_at_time(reader, result, top_node_nbrs, -1e-4, 6, time=1.0)

    def _assert_oedometer_effective_stresses(self, effective_stresses, expected_yy, places_yy):
        for element in effective_stresses:
            for integration_point in element:
                self.assertAlmostEqual(0.0, integration_point[0,0], 3)
                self.assertAlmostEqual(expected_yy, integration_point[1,1], places_yy)
                self.assertAlmostEqual(0.0, integration_point[2,2], 3)

    def _assert_y_displacements_at_time(self, reader, result, node_ids, expected_y, places, time=1.0):
        displacements = reader.nodal_values_at_time("DISPLACEMENT", time, result, node_ids=node_ids)
        for displacement in displacements:
            self.assertAlmostEqual(expected_y, displacement[1], places)

    def _assert_integration_point_tensor_results(self, integration_point_tensors, expected_integration_point_tensor, places, result_name):
        for idx, ip_tensor in enumerate(integration_point_tensors):
            self.assertVectorAlmostEqual(expected_integration_point_tensor, ip_tensor, places, msg = f"{result_name} components at integration point {idx}")

    def _assert_integration_point_tensors(self, reader, result, variable_name, expected_tensors, places, time):
        """Assert tensor values for all integration points across both elements."""
        for flat_index, expected_tensor in enumerate(expected_tensors):
            element_id, ip_index = self._get_element_id_and_ip_index(flat_index)
            tensor = reader.element_integration_point_values_at_time(variable_name, time, result, [element_id], [ip_index])[0]
            self._assert_integration_point_tensor_results(tensor, expected_tensor, places, variable_name)

    def _get_element_id_and_ip_index(self, flat_index):
        """Map a flat integration-point index to the corresponding element and local integration point index."""
        if flat_index < 3:
            return 1, flat_index
        return 2, flat_index - 3

    def _assert_first_ip_tensor(self, reader, result, variable_name, expected_values, precision_places, time=1.0):
        for element_id, expected_tensor in enumerate(expected_values, start=1):
            tensor = reader.element_integration_point_values_at_time(variable_name, time, result, [element_id], [0])[0]
            self._assert_integration_point_tensor_results(tensor, expected_tensor, precision_places, variable_name)

if __name__ == '__main__':
    KratosUnittest.main()
