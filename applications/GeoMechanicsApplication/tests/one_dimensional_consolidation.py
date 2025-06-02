import math
import os
import shutil

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest

import test_helper
import analytical_solutions


class OneDimensionalConsolidationTestBase(KratosUnittest.TestCase):
    def setUp(self):
        super().setUp()

        self.test_root = test_helper.get_file_path("one_dimensional_consolidation")
        self.test_path = os.path.join(self.test_root, self._get_test_dir_name())

        shutil.rmtree(self.test_path, ignore_errors=True)

        os.makedirs(self.test_path)

        self.number_of_stages = 11
        self.project_parameters_filenames = [f"ProjectParameters_stage{i+1}.json" for i in range(self.number_of_stages)]
        input_filenames = self.project_parameters_filenames + ["MaterialParameters.json", "1D-Consolidationtest.mdpa"]

        for filename in input_filenames:
            shutil.copy(os.path.join(self.test_root, filename), os.path.join(self.test_path, filename))

        self.top_node_ids = [197, 198, 199, 200, 201]

        self.end_times = [8640, 17280, 43200, 86400, 172800, 432000, 864000, 1728000, 4320000, 8640000]

        c_v = 1.0 / 864000.0  # consolidation coefficient
        self.h = 1.0  # sample height
        self.t_vs = [c_v * t / (self.h * self.h) for t in self.end_times]


    def _get_test_dir_name(self):
        raise RuntimeError("This base class does not provide a generic test directory name")


    def _get_y_coordinates_of_all_nodes(self):
        post_msh_file_path = os.path.join(self.test_path, "1D-Consolidationtest_stage1.post.msh")
        return [coord[1] + 1.0 for coord in test_helper.read_coordinates_from_post_msh_file(post_msh_file_path)]


    def _get_analytical_relative_water_pressures(self, t_v, y_coordinates):
        return [analytical_solutions.calculate_relative_water_pressure(y, self.h, t_v) for y in y_coordinates]


    def _get_numerical_relative_water_pressures(self, time, stage_no):
        output_file_path = os.path.join(self.test_path, f"1D-Consolidationtest_stage{stage_no}.post.res")
        reader = test_helper.GiDOutputFileReader()
        output_data = reader.read_output_from(output_file_path)
        # Invert the sign of the water pressures resulting from the numerical solution, to make them match the
        # analytical solution which assumes compressive water pressures to be positive rather than negative.
        # In addition, since the load has size 1, the relative values are equal to the total values.
        return [-1.0 * pw for pw in reader.nodal_values_at_time("WATER_PRESSURE", time, output_data)]


    def _get_numerical_settlement_values(self, time, stage_no, node_ids):
        output_file_path = os.path.join(self.test_path, f"1D-Consolidationtest_stage{stage_no}.post.res")
        reader = test_helper.GiDOutputFileReader()
        output_data = reader.read_output_from(output_file_path)
        # Settlement corresponds to downward total displacement. Invert the sign of the total displacement values to
        # make them match the analytical solution which assumes the settlement values to be positive.
        return [-1.0 * u[1] for u in reader.nodal_values_at_time("TOTAL_DISPLACEMENT", time, output_data, node_ids)]


    def _calculate_rmse_of_differences(self, values1, values2):
        differences = [value1 - value2 for value1, value2 in zip(values1, values2)]
        return math.sqrt(sum([diff * diff for diff in differences]) / len(differences))


    def _check_relative_water_pressures(self):
        y_coordinates = self._get_y_coordinates_of_all_nodes()
        stage_no = 2  # The first stage is not checked
        rmse_values = []
        for t_v, time in zip(self.t_vs, self.end_times):
            analytical_solution = self._get_analytical_relative_water_pressures(t_v, y_coordinates)
            numerical_solution = self._get_numerical_relative_water_pressures(time, stage_no)
            rmse_values.append(self._calculate_rmse_of_differences(numerical_solution, analytical_solution))
            stage_no += 1

        self._check_rmse_values(rmse_values, "relative water pressure values")


    def _check_degree_of_consolidation(self):
        # calculate the degree of consolidation analytical solution for all stages and calculate the error
        # Verruijt's notations
        q = -1.0                            # the load
        K = 3.33e2                          # the compression modulus
        G = 5.0e2                           # shear modulus
        m_v = 1/(K + 4/3 * G)               # the compressibility coefficient
        beta = 0.5e-9                       # the compressibility of the water.
        n = 0.3                             # the porosity
        delta_h0 = (-1) * m_v * self.h * q * n * beta / (m_v + n * beta) # deformation immediately after the application of the load
        delta_h_infinity = (-1) * m_v * self.h * q # the final deformation

        stage_no = 2  # The first stage is not checked
        rmse_values = []
        for t_v, time in zip(self.t_vs, self.end_times):
            settlement_values = self._get_numerical_settlement_values(time, stage_no, self.top_node_ids)
            numerical_degree_values = [(u_y - delta_h0) / (delta_h_infinity - delta_h0) for u_y in settlement_values]
            analytical_degree = analytical_solutions.calculate_degree_of_1d_consolidation(t_v)
            rmse_values.append(self._calculate_rmse_of_differences(numerical_degree_values, [analytical_degree] * len(self.top_node_ids)))
            stage_no += 1

        self._check_rmse_values(rmse_values, "degree of consolidation values")


    def _check_rmse_values(self, values, description):
        accuracy = 0.01
        for stage_index, value in enumerate(values):
            self.assertLess(value, accuracy, msg=f"RMSE of {description} in stage {stage_index+2}")  # Checking starts from stage 2

class KratosGeoMechanics1DConsolidation(OneDimensionalConsolidationTestBase):
    """
    This class contains benchmark tests which are checked with the analytical solution
    """
    def _get_test_dir_name(self):
        return "python"


    def test_1d_consolidation(self):
        test_helper.run_stages(self.test_path, self.number_of_stages)

        self._check_relative_water_pressures()
        self._check_degree_of_consolidation()


class KratosGeoMechanics1DConsolidationCppRoute(OneDimensionalConsolidationTestBase):
    def _get_test_dir_name(self):
        return "cpp"


    def test_1d_consolidation(self):
        import KratosMultiphysics.GeoMechanicsApplication.run_geo_settlement as run_geo_settlement

        status = run_geo_settlement.run_stages(self.test_path, self.project_parameters_filenames)
        self.assertEqual(status, 0)

        self._check_relative_water_pressures()
        self._check_degree_of_consolidation()


if __name__ == '__main__':
    KratosUnittest.main()
