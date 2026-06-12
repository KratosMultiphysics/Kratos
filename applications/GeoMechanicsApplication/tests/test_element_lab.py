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
        file_path = test_helper.get_file_path(os.path.join('test_element_lab', 'test_triaxial', 'drained'))
        expected_disp = test_helper.get_values_from_csv_as_vectors(
            Path(file_path) / "expected_disp.csv",
            ["node_id"],
            ["disp_x", "disp_y", "disp_z"])
        expected_stress = test_helper.get_values_from_csv_as_vectors(
            Path(file_path) / "expected_stress.csv",
            ["element_id", "ip_index"],
            ["stress_xx", "stress_yy", "stress_zz", "stress_xy", "stress_yz", "stress_xz"])
        expected_strain = test_helper.get_values_from_csv_as_vectors(
            Path(file_path) / "expected_strain.csv",
            ["element_id", "ip_index"],
            ["strain_xx", "strain_yy", "strain_zz", "strain_xy", "strain_yz", "strain_xz"])
        self._run_triaxial_regression_test(file_path, 'triaxial_test_output.post.res', expected_disp, expected_stress,
                                           expected_strain, 4)

    def test_triaxial_undrained(self):
        """
        Regression test for the undrained triaxial experiment.
        """
        file_path = test_helper.get_file_path(os.path.join('test_element_lab', 'test_triaxial', 'undrained'))
        expected_disp = test_helper.get_values_from_csv_as_vectors(
            Path(file_path) / "expected_disp.csv",
            ["node_id"],
            ["disp_x", "disp_y", "disp_z"])
        expected_stress = test_helper.get_values_from_csv_as_vectors(
            Path(file_path) / "expected_stress.csv",
            ["element_id", "ip_index"],
            ["stress_xx", "stress_yy", "stress_zz", "stress_xy", "stress_yz", "stress_xz"])
        expected_strain = test_helper.get_values_from_csv_as_vectors(
            Path(file_path) / "expected_strain.csv",
            ["element_id", "ip_index"],
            ["strain_xx", "strain_yy", "strain_zz", "strain_xy", "strain_yz", "strain_xz"])
        self._run_triaxial_regression_test(file_path, 'triaxial_undrained_test_output.post.res', expected_disp,
                                           expected_stress, expected_strain, 4)

    def _run_triaxial_regression_test(self, file_path, output_file_name, expected_displacement, expected_stress, expected_strain, precision_places):
        test_helper.run_kratos(file_path)
        reader = GiDOutputFileReader()
        result = reader.read_output_from(os.path.join(file_path, output_file_name))
        time = 1.0
        node_ids = list(expected_displacement.keys())
        expected_disp = [expected_displacement[node_id] for node_id in node_ids]
        self.assert_nodal_values_at_time(result, "DISPLACEMENT", node_ids, expected_disp, time, precision_places)
        self.assert_integration_point_tensors(result, "CAUCHY_STRESS_TENSOR", expected_stress, time, precision_places)
        self.assert_integration_point_tensors(result, "ENGINEERING_STRAIN_TENSOR", expected_strain, time, precision_places)
            
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
        expected_water_pressures = [[0.0] * 9 for _ in range(nr_of_phases)]

        linear_elastic_file_path = test_helper.get_file_path(os.path.join('test_element_lab', 'test_crs', stage_name, 'linear_elastic'))
        linear_elastic_stresses = test_helper.get_values_from_csv_grouped(
            Path(linear_elastic_file_path) / "expected_stress.csv",
            ['phase'], ['element_id', 'ip_index'],
            ["stress_xx", "stress_yy", "stress_zz", "stress_yz", "stress_xz", "stress_xy"])
        linear_elastic_strains = test_helper.get_values_from_csv_grouped(
            Path(linear_elastic_file_path) / "expected_strain.csv",
            ['phase'], ['element_id', 'ip_index'],
            ["strain_xx", "strain_yy", "strain_zz", "strain_yz", "strain_xz", "strain_xy"])
        expected_stresses = [linear_elastic_stresses[phase] for phase in range(1, nr_of_phases + 1)]
        expected_strains = [linear_elastic_strains[phase] for phase in range(1, nr_of_phases + 1)]

        self._run_crs_test_for_model(stage_name, "linear_elastic")
        self._check_crs_results(stage_name, "linear_elastic", nr_of_phases, expected_stresses, expected_strains, expected_water_pressures, delta_stress=1.0)

        mohr_coulomb_file_path = test_helper.get_file_path(os.path.join('test_element_lab', 'test_crs', stage_name, 'mohr_coulomb'))
        mohr_coulomb_stresses = test_helper.get_values_from_csv_grouped(
            Path(mohr_coulomb_file_path) / "expected_stress.csv",
            ['phase'], ['element_id', 'ip_index'],
            ["stress_xx", "stress_yy", "stress_zz", "stress_xy", "stress_yz", "stress_xz"])
        mohr_coulomb_strains = test_helper.get_values_from_csv_grouped(
            Path(mohr_coulomb_file_path) / "expected_strain.csv",
            ['phase'], ['element_id', 'ip_index'],
            ["strain_xx", "strain_yy", "strain_zz", "strain_xy", "strain_yz", "strain_xz"])
        expected_stresses = [mohr_coulomb_stresses[phase] for phase in range(1, nr_of_phases + 1)]
        expected_strains = [mohr_coulomb_strains[phase] for phase in range(1, nr_of_phases + 1)]

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
        node_ids = list(range(1, 10))
        for i in range(nr_of_phases):
            self.assert_integration_point_tensors(
                result, "CAUCHY_STRESS_TENSOR", expected_stresses[i], times[i], places=places_stress, delta=delta_stress)
            self.assert_integration_point_tensors(
                result, "ENGINEERING_STRAIN_TENSOR", expected_strains[i], times[i], places=places_strain)
            self.assert_nodal_values_at_time(result, "WATER_PRESSURE", node_ids, expected_water_pressures[i], times[i], places=places_water_pressure)

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
