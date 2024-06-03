import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper


class KratosGeoMechanicsCompressibilityTests(KratosUnittest.TestCase):
    """
    This class
    """

    def test_compressibility_upw_small_strain(self):
        test_name = 'compressibility_tests/upw_small_strain'
        file_path = test_helper.get_file_path(test_name)
        test_helper.run_kratos(file_path)

        output_file_path = os.path.join(file_path, 'output.post.res')
        output_reader = test_helper.GiDOutputFileReader()
        output_data = output_reader.read_output_from(output_file_path)

        top_node_nbrs = [5, 6, 7, 8]
        water_pressures_at_top = test_helper.GiDOutputFileReader.nodal_values_at_time("WATER_PRESSURE", 1.0, output_data,
                                                                                    node_ids=top_node_nbrs)
        bottom_node_nbrs = [1, 2, 3, 4]
        water_pressures_at_bottom = test_helper.GiDOutputFileReader.nodal_values_at_time("WATER_PRESSURE", 1.0, output_data,
                                                                                       node_ids=bottom_node_nbrs)

        for water_pressure_at_top in water_pressures_at_top:
            self.assertAlmostEqual(water_pressure_at_top, 0.0, 6)
        for water_pressure_at_bottom in water_pressures_at_bottom:
            self.assertAlmostEqual(water_pressure_at_bottom, -49.1786, 6)

    def test_compressibility_upw_small_strain_fic(self):
        test_name = 'compressibility_tests/upw_small_strain_fic'
        file_path = test_helper.get_file_path(test_name)
        test_helper.run_kratos(file_path)

        output_file_path = os.path.join(file_path, 'output.post.res')
        output_reader = test_helper.GiDOutputFileReader()
        output_data = output_reader.read_output_from(output_file_path)

        top_node_nbrs = [5, 6, 7, 8]
        water_pressures_at_top = test_helper.GiDOutputFileReader.nodal_values_at_time("WATER_PRESSURE", 1.0, output_data,
                                                                                      node_ids=top_node_nbrs)
        bottom_node_nbrs = [1, 2, 3, 4]
        water_pressures_at_bottom = test_helper.GiDOutputFileReader.nodal_values_at_time("WATER_PRESSURE", 1.0, output_data,
                                                                                         node_ids=bottom_node_nbrs)

        for water_pressure_at_top in water_pressures_at_top:
            self.assertAlmostEqual(water_pressure_at_top, 0.0, 6)
        for water_pressure_at_bottom in water_pressures_at_bottom:
            self.assertAlmostEqual(water_pressure_at_bottom, -57.693, 6)

    def test_compressibility_upw_small_strain_updated_lagrange_fic(self):
        test_name = 'compressibility_tests/upw_small_strain_updated_lagrange_fic'
        file_path = test_helper.get_file_path(test_name)
        test_helper.run_kratos(file_path)

        output_file_path = os.path.join(file_path, 'output.post.res')
        output_reader = test_helper.GiDOutputFileReader()
        output_data = output_reader.read_output_from(output_file_path)

        top_node_nbrs = [5, 6, 7, 8]
        water_pressures_at_top = test_helper.GiDOutputFileReader.nodal_values_at_time("WATER_PRESSURE", 1.0, output_data,
                                                                                      node_ids=top_node_nbrs)
        bottom_node_nbrs = [1, 2, 3, 4]
        water_pressures_at_bottom = test_helper.GiDOutputFileReader.nodal_values_at_time("WATER_PRESSURE", 1.0, output_data,
                                                                                         node_ids=bottom_node_nbrs)

        for water_pressure_at_top in water_pressures_at_top:
            self.assertAlmostEqual(water_pressure_at_top, 0.0, 6)
        for water_pressure_at_bottom in water_pressures_at_bottom:
            self.assertAlmostEqual(water_pressure_at_bottom, -57.5806, 6)


    def test_compressibility_upw_small_strain_updated_lagrange(self):
        test_name = 'compressibility_tests/upw_small_strain_updated_lagrange'
        file_path = test_helper.get_file_path(test_name)
        test_helper.run_kratos(file_path)

        output_file_path = os.path.join(file_path, 'output.post.res')
        output_reader = test_helper.GiDOutputFileReader()
        output_data = output_reader.read_output_from(output_file_path)

        top_node_nbrs = [5, 6, 7, 8]
        water_pressures_at_top = test_helper.GiDOutputFileReader.nodal_values_at_time("WATER_PRESSURE", 1.0, output_data,
                                                                                      node_ids=top_node_nbrs)
        bottom_node_nbrs = [1, 2, 3, 4]
        water_pressures_at_bottom = test_helper.GiDOutputFileReader.nodal_values_at_time("WATER_PRESSURE", 1.0, output_data,
                                                                                         node_ids=bottom_node_nbrs)

        for water_pressure_at_top in water_pressures_at_top:
            self.assertAlmostEqual(water_pressure_at_top, 0.0, 6)
        for water_pressure_at_bottom in water_pressures_at_bottom:
            self.assertAlmostEqual(water_pressure_at_bottom, -49.0716, 6)


if __name__ == '__main__':
    KratosUnittest.main()
