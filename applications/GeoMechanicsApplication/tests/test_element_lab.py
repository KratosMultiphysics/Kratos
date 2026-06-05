import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import GiDOutputFileReader
import KratosMultiphysics.GeoMechanicsApplication.run_multiple_stages as run_multiple_stages
import test_helper
import KratosGeoUnittest
from pathlib import Path

class KratosGeoMechanicsLabElementTests(KratosGeoUnittest.TestCase):
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
            self.assert_integration_point_tensors(result, "CAUCHY_STRESS_TENSOR",
                self._make_integration_point_tensor_entries(expected_stress, num_elements=2, num_integration_points_per_element=3),
                time ,precision_places)
            self.assert_integration_point_tensors(result, "ENGINEERING_STRAIN_TENSOR",
                self._make_integration_point_tensor_entries(expected_strain, num_elements=2, num_integration_points_per_element=3),
                time, precision_places)
        else:
            self.assert_integration_point_tensors(result, "CAUCHY_STRESS_TENSOR", 
                self._make_integration_point_tensor_entries(expected_stress, num_elements=2, num_integration_points_per_element=1),
                time, precision_places)
            self.assert_integration_point_tensors(result, "ENGINEERING_STRAIN_TENSOR", 
                self._make_integration_point_tensor_entries(expected_strain, num_elements=2, num_integration_points_per_element=1),
                time, precision_places)
            
    def test_oedometer_drained(self):
        """Regression test for the oedometer experiment on a linear elastic model with constant pore water pressure."""
        expected_stress_per_ip = [-1e+06/3.0, -1e+06, -1e+06/3.0, 0.0, 0.0, 0.0]
        expected_stress = [expected_stress_per_ip] * 6
        expected_y_displacements = [-0.0208, -0.0417, -0.0625, -0.0833]
        displacement_times = [0.25, 0.5, 0.75, 1.0]
        top_nodes = [7, 8, 9]
        self._run_oedometer_regression_test('drained', 'test_oedometer_output.post.res', expected_stress,
                                            expected_y_displacements, displacement_times, top_nodes, 0, 4)
        
    def _run_oedometer_regression_test(self, stage_name, output_file_name, expected_stress,
                                       expected_y_displacements, displacement_times, top_nodes, precision_places_stress, precision_places_displacement=4):
        """Run the oedometer regression test and validate stress tensors and y-displacements in time."""
        test_name = 'test_oedometer'
        file_path = test_helper.get_file_path(os.path.join('test_element_lab', test_name, stage_name))
        test_helper.run_kratos(file_path)

        reader = GiDOutputFileReader()
        result = reader.read_output_from(os.path.join(file_path, output_file_name))

        self.assert_integration_point_tensors(
            result, "CAUCHY_STRESS_TENSOR",
            self._make_integration_point_tensor_entries(expected_stress, num_elements=2, num_integration_points_per_element=3),
            1.0, precision_places_stress)

        for time, expected_y in zip(displacement_times, expected_y_displacements):
            self.assert_y_displacements_at_time(result, top_nodes, expected_y, time, precision_places_displacement)

    def test_dss_drained(self):
        """Regression test for the direct simple shear experiment with constant pore water pressure."""
        stage_name = 'drained'
        expected_stress = [[-1.0e+05, -1.0e+05, -1.0e+05, 8.0e+05, 0.0, 0.0]] * 6
        expected_strain = [[0.0, 0.0, 0.0, 0.1, 0.0, 0.0]] * 6
        self._run_dss_regression_test(stage_name, "linear_elastic", "test_dss_output.post.res", expected_stress, 3, expected_strain, 6, time = 1.0)

        expected_stress = [[-205490, -205490, -152745, 88656.5, 0.0, 0.0]] * 6
        expected_strain = [[0.0, 0.0, 0.0, 0.1, 0.0, 0.0]] * 6
        self._run_dss_regression_test(stage_name, "mohr_coulomb", "test_dss_output.post.res", expected_stress, 3, expected_strain, 6, time = 1.0)
    
    def _run_dss_regression_test(self, stage_name, model_name, output_file_name, expected_stress, places_stress, expected_strain, places_strain, time=1.0):
        """Run the direct simple shear regression test and validate stress and strain tensors at given times."""
        test_name = 'test_dss'
        common_file_path = test_helper.get_file_path(os.path.join('test_element_lab', test_name, "common"))
        material_file_path = test_helper.get_file_path(os.path.join('test_element_lab', test_name, stage_name, model_name))
        run_multiple_stages.run_stages(material_file_path, 1, filename_pattern="ProjectParameters.json", input_path=common_file_path)

        reader = GiDOutputFileReader()
        result = reader.read_output_from(os.path.join(material_file_path, output_file_name))

        self.assert_integration_point_tensors(
            result, "CAUCHY_STRESS_TENSOR",
            self._make_integration_point_tensor_entries(expected_stress, num_elements=2, num_integration_points_per_element=3),
            time, places_stress)
        self.assert_integration_point_tensors(
            result, "ENGINEERING_STRAIN_TENSOR",
            self._make_integration_point_tensor_entries(expected_strain, num_elements=2, num_integration_points_per_element=3),
            time, places_strain)
        
    def test_crs_drained(self):
        """Regression test for the CRS experiment with constant pore water pressure."""
        stage_name = 'drained'
        nr_of_phases = 5
        nr_of_elements = 2
        nr_of_integration_points_per_element = 3
        total_integration_points = nr_of_elements * nr_of_integration_points_per_element

        expected_strains = [
            self._repeat_tensor([0.0, -0.1, 0.0, 0.0, 0.0, 0.0], total_integration_points),
            self._repeat_tensor([0.0, -0.05, 0.0, 0.0, 0.0, 0.0], total_integration_points),
            self._repeat_tensor([0.0, -0.25, 0.0, 0.0, 0.0, 0.0], total_integration_points),
            self._repeat_tensor([0.0, -0.25, 0.0, 0.0, 0.0, 0.0], total_integration_points),
            self._repeat_tensor([0.0, -0.40, 0.0, 0.0, 0.0, 0.0], total_integration_points),
        ]
        expected_water_pressures = [[0.0] * 9 for _ in range(nr_of_phases)]

        expected_stresses = [
            self._repeat_tensor([-4e+05, -12e+05, -4e+05, 0.0, 0.0, 0.0], total_integration_points),
            self._repeat_tensor([-2e+05, -6e+05, -2e+05, 0.0, 0.0, 0.0], total_integration_points),
            self._repeat_tensor([-1e+06, -3e+06, -1e+06, 0.0, 0.0, 0.0], total_integration_points),
            self._repeat_tensor([-1e+06, -3e+06, -1e+06, 0.0, 0.0, 0.0], total_integration_points),
            self._repeat_tensor([-1.6e+06, -4.8e+06, -1.6e+06, 0.0, 0.0, 0.0], total_integration_points),
        ]

        self._run_crs_test_for_model(stage_name, "linear_elastic")
        self._check_crs_results(stage_name, "linear_elastic", nr_of_phases, expected_stresses, expected_strains, expected_water_pressures, delta_stress=1.0)

        # mohr coulomb expectations per phase
        expected_stresses = [
            self._repeat_tensor([-9.01878e+05, -2.22843e+06, -9.01878e+05, 0.0, 0.0, 0.0], total_integration_points),
            self._repeat_tensor([-5.01878e+05, -1.02843e+06, -5.01878e+05, 0.0, 0.0, 0.0], total_integration_points),
            self._repeat_tensor([-2.25697e+06, -5.56725e+06, -2.25697e+06, 0.0, 0.0, 0.0], total_integration_points),
            self._repeat_tensor([-2.25697e+06, -5.56725e+06, -2.25697e+06, 0.0, 0.0, 0.0], total_integration_points),
            self._repeat_tensor([-3.61205e+06, -8.90607e+06, -3.61205e+06, 0.0, 0.0, 0.0], total_integration_points)
            ]

        self._run_crs_test_for_model(stage_name, "mohr_coulomb")
        self._check_crs_results(stage_name, "mohr_coulomb", nr_of_phases, expected_stresses, expected_strains, expected_water_pressures, delta_stress=1.0)

    def _run_crs_test_for_model(self, stage_name, model_name):
        test_name =  "test_crs"
        common_file_path = test_helper.get_file_path(os.path.join('test_element_lab', test_name, "common"))
        material_file_path = test_helper.get_file_path(os.path.join('test_element_lab', test_name, stage_name, model_name))
        run_multiple_stages.run_stages(material_file_path, 1, filename_pattern="ProjectParameters.json", input_path=common_file_path)
    
    def _check_crs_results(self, stage_name, model_name, nr_of_phases, expected_stresses, expected_strains, expected_water_pressures, times = [3600.0, 7200.0, 10800.0, 14400.0, 18000.0], places_stress=None, places_strain=4, places_water_pressure=2, delta_stress = None):
        reader = GiDOutputFileReader()
        file_path = test_helper.get_file_path(Path('test_element_lab') / "test_crs" / stage_name / model_name)
        result = reader.read_output_from(Path(file_path) / "test_crs_output.post.res")
        for i in range(nr_of_phases):
            self.assert_integration_point_tensors(
                result, "CAUCHY_STRESS_TENSOR",
                self._make_integration_point_tensor_entries(expected_stresses[i], num_elements=2, num_integration_points_per_element=3),
                times[i], places=places_stress, delta=delta_stress)
            self.assert_integration_point_tensors(
                result, "ENGINEERING_STRAIN_TENSOR",
                self._make_integration_point_tensor_entries(expected_strains[i], num_elements=2, num_integration_points_per_element=3),
                times[i], places=places_strain)
            self.assert_nodal_values_at_time(result, "WATER_PRESSURE", expected_water_pressures[i], times[i], places=places_water_pressure)

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
        self.assert_integration_point_tensors(
            result, "CAUCHY_STRESS_TENSOR",
            self._make_integration_point_tensor_entries(expected_stress, num_elements=2, num_integration_points_per_element=3),
            time=1.0, delta=0.002)

        result = reader.read_output_from(os.path.join(project_path, "triaxial_comp_6n_stage2.post.res"))
        expected_stress = [[-100.0, -300.0, -100.0, 0.0, 0.0, 0.0]] * 6
        self.assert_integration_point_tensors(
            result, "CAUCHY_STRESS_TENSOR",
            self._make_integration_point_tensor_entries(expected_stress, num_elements=2, num_integration_points_per_element=3),
            time=1.25, delta=0.002)

    def test_oedometer_ULFEM(self):
        """
        Oedometer test on a linear elastic model with 2D6N elements
        """
        test_name = 'oedometer_ULFEM'
        file_path = test_helper.get_file_path(os.path.join('test_element_lab', test_name))
        simulation = test_helper.run_kratos(file_path)
        effective_stresses = test_helper.get_cauchy_stress_tensor(simulation)
        self._assert_oedometer_effective_stresses(effective_stresses, -1e+06, 2)

        top_node_nbrs = [1, 2]
        output_file = os.path.join(file_path, test_name+'.post.res')
        reader = GiDOutputFileReader()
        result = reader.read_output_from(output_file)

        self.assert_y_displacements_at_time(result, top_node_nbrs, -0.0099, 0.1, 4)
        self.assert_y_displacements_at_time(result, top_node_nbrs, -0.0654, 0.7, 4)
        self.assert_y_displacements_at_time(result, top_node_nbrs, -0.0909, 1.0, 4)

    def test_oedometer_ULFEM_diff_order(self):
        """
        Oedometer test on a linear elastic model with 2D6N with different order elements
        """
        test_name = 'oedometer_ULFEM_diff_order'
        project_path = test_helper.get_file_path(os.path.join('test_element_lab', test_name))
        simulation = test_helper.run_kratos(project_path)
        effective_stresses = test_helper.get_cauchy_stress_tensor(simulation)
        self._assert_oedometer_effective_stresses(effective_stresses, -1e+03, 3)

        output_file = os.path.join(project_path, test_name + '.post.res')
        reader = GiDOutputFileReader()
        result = reader.read_output_from(output_file)
        top_node_nbrs = [1]
        self.assert_y_displacements_at_time(result, top_node_nbrs, -1e-4, 1.0, 6)

    def _assert_oedometer_effective_stresses(self, effective_stresses, expected_yy, places_yy):
        for element in effective_stresses:
            for integration_point in element:
                self.assertAlmostEqual(0.0, integration_point[0,0], 3)
                self.assertAlmostEqual(expected_yy, integration_point[1,1], places_yy)
                self.assertAlmostEqual(0.0, integration_point[2,2], 3)

    def _repeat_tensor(self, tensor, count):
        return [list(tensor) for _ in range(count)]

    def _make_integration_point_tensor_entries(self, expected_tensors, num_elements, num_integration_points_per_element):
        result = {}
        idx = 0
        for element_id in range(1, num_elements + 1):
            for ip_index in range(num_integration_points_per_element):
                result[(element_id, ip_index)] = expected_tensors[idx]
                idx += 1
        return result 

if __name__ == '__main__':
    KratosUnittest.main()
