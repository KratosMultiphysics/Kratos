import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper


class KratosGeoMechanicsCompressibilityTests(KratosUnittest.TestCase):
    """
    This class contains simple single-element regression tests, checking the water pressure after a single time step.
    It is a workflow test, which will fail if the additions of right hand side contributions (e.g. compressibility,
    permeability, coupling) are changed.
    """

    def test_compressibility_upw_small_strain(self):
        test_name = 'compressibility_tests/upw_small_strain'
        expected_water_pressure_at_bottom = -49.2298
        self.run_and_assert_water_pressures(expected_water_pressure_at_bottom, test_name)

    def test_compressibility_upw_small_strain_fic(self):
        test_name = 'compressibility_tests/upw_small_strain_fic'
        expected_water_pressure_at_bottom = -56.1231
        self.run_and_assert_water_pressures(expected_water_pressure_at_bottom, test_name)

    def test_compressibility_upw_updated_lagrange_fic(self):
        test_name = 'compressibility_tests/upw_updated_lagrange_fic'
        expected_water_pressure_at_bottom = -55.6957
        self.run_and_assert_water_pressures(expected_water_pressure_at_bottom, test_name)

    def test_compressibility_upw_updated_lagrange(self):
        test_name = 'compressibility_tests/upw_updated_lagrange'
        expected_water_pressure_at_bottom = -48.8369
        self.run_and_assert_water_pressures(expected_water_pressure_at_bottom, test_name)

    def test_compressibility_upw_small_strain_diff_order(self):
        test_name = 'compressibility_tests/upw_small_strain_diff_order'
        expected_water_pressures_at_bottom = [-49.4291, -49.4291, -49.4291, -49.4291]
        self.run_and_assert_water_pressures_diff_order(expected_water_pressures_at_bottom, test_name)

    def test_compressibility_upw_diff_order_updated_lagrange(self):
        test_name = 'compressibility_tests/upw_diff_order_updated_lagrange'
        expected_water_pressures_at_bottom = [-48.8613, -48.8613, -48.8613, -48.8613]
        self.run_and_assert_water_pressures_diff_order(expected_water_pressures_at_bottom, test_name)

    def run_and_assert_water_pressures(self, expected_water_pressure_at_bottom, test_name):
        output_data = self.run_simulation(test_name)
        water_pressures_at_bottom, water_pressures_at_top = self.get_water_pressures(output_data)
        for water_pressure_at_top in water_pressures_at_top:
            self.assertAlmostEqual(water_pressure_at_top, 0.0, 6)
        for water_pressure_at_bottom in water_pressures_at_bottom:
            self.assertAlmostEqual(water_pressure_at_bottom, expected_water_pressure_at_bottom, 4)

    def run_and_assert_water_pressures_diff_order(self, expected_water_pressures_at_bottom, test_name):
        output_data = self.run_simulation(test_name)

        water_pressures_at_bottom, water_pressures_at_top = self.get_water_pressures(output_data)
        for water_pressure_at_top in water_pressures_at_top:
            self.assertAlmostEqual(water_pressure_at_top, 0.0, 6)
        for water_pressure_at_bottom, expected_water_pressure_at_bottom in (
                zip(water_pressures_at_bottom, expected_water_pressures_at_bottom)):
            self.assertAlmostEqual(water_pressure_at_bottom, expected_water_pressure_at_bottom, 4)

    @staticmethod
    def run_simulation(test_name):
        file_path = test_helper.get_file_path(test_name)
        test_helper.run_kratos(file_path)
        output_file_path = os.path.join(file_path, 'output.post.res')
        output_reader = test_helper.GiDOutputFileReader()
        return output_reader.read_output_from(output_file_path)

    @staticmethod
    def get_water_pressures(output_data):
        top_node_nbrs = [5, 6, 7, 8]
        water_pressures_at_top = test_helper.GiDOutputFileReader.nodal_values_at_time("WATER_PRESSURE", 1.0,
                                                                                      output_data,
                                                                                      node_ids=top_node_nbrs)
        bottom_node_nbrs = [1, 2, 3, 4]
        water_pressures_at_bottom = test_helper.GiDOutputFileReader.nodal_values_at_time("WATER_PRESSURE", 1.0,
                                                                                         output_data,
                                                                                         node_ids=bottom_node_nbrs)
        return water_pressures_at_bottom, water_pressures_at_top


if __name__ == '__main__':
    KratosUnittest.main()
