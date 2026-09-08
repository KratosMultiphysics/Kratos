import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import GiDOutputFileReader
import test_helper


class KratosGeoMechanicsFixedSpatialVariationTests(KratosUnittest.TestCase):

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def test_fixed_spatial_variation(self):
        test_name = 'fixed_spatial_variation'
        file_path = test_helper.get_file_path(test_name)
        test_helper.run_kratos(file_path)
        result_file_name = os.path.join(file_path, test_name + '.post.res')
        reader = GiDOutputFileReader()
        actual_data = reader.read_output_from(result_file_name)

        top_nodes = [3, 4, 7]
        actual_top_nodal_pressures = reader.nodal_values_at_time("WATER_PRESSURE", 1.0, actual_data, top_nodes)
        expected_top_pressures = [-5000, 0.0, -2500]
        for actual_pressure, expected_pressure in zip(actual_top_nodal_pressures, expected_top_pressures):
            self.assertAlmostEqual(actual_pressure, expected_pressure, 3)

        bottom_nodes = [1, 2, 5]
        actual_bottom_nodal_pressures = reader.nodal_values_at_time("WATER_PRESSURE", 1.0, actual_data, bottom_nodes)
        expected_bottom_pressures = [-40000, -30000, -35000]
        for actual_pressure, expected_pressure in zip(actual_bottom_nodal_pressures, expected_bottom_pressures):
            self.assertAlmostEqual(actual_pressure, expected_pressure, 3)


if __name__ == '__main__':
    KratosUnittest.main()
