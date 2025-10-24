import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper


class KratosGeoMechanicsTransientThermalValidationTests(KratosUnittest.TestCase):
    """
    This class contains benchmark tests which are checked with the regression on a previously obtained value.
    """
    etalon_value1 = 28.04411163544510063559
    etalon_value2 = 17.55892791313559322
    etalon_value3 = 41.3797035928672316
    etalon_value4 = 35.31073446327683615819
    etalon_value5 = 26.6666666666666667
    etalon_value6 = 88.2768361581920904

    def simulate_thermal_case(self, test_name):
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', test_name))
        simulation = test_helper.run_kratos(file_path)
        return test_helper.get_temperature(simulation)

    def test_thermal_heat_flux_2D6N(self):
        temperature = self.simulate_thermal_case('test_thermal_heat_flux/test_thermal_heat_flux_2D6N')
        self.assertAlmostEqual(self.etalon_value1, temperature[57])

    def test_thermal_heat_flux_2D10N(self):
        temperature = self.simulate_thermal_case('test_thermal_heat_flux/test_thermal_heat_flux_2D10N')
        self.assertAlmostEqual(self.etalon_value1, temperature[77])

    def test_thermal_heat_flux_2D15N(self):
        temperature = self.simulate_thermal_case('test_thermal_heat_flux/test_thermal_heat_flux_2D15N')
        self.assertAlmostEqual(self.etalon_value1, temperature[97])

    def test_thermal_heat_flux_2D4N(self):
        temperature = self.simulate_thermal_case('test_thermal_heat_flux/test_thermal_heat_flux_2D4N')
        self.assertAlmostEqual(self.etalon_value2, temperature[18])

    def test_thermal_heat_flux_2D8N(self):
        temperature = self.simulate_thermal_case('test_thermal_heat_flux/test_thermal_heat_flux_2D8N')
        self.assertAlmostEqual(self.etalon_value2, temperature[50])

    def test_thermal_heat_flux_2D9N(self):
        temperature = self.simulate_thermal_case('test_thermal_heat_flux/test_thermal_heat_flux_2D9N')
        self.assertAlmostEqual(self.etalon_value2, temperature[63])

    def test_thermal_heat_flux_3D10N(self):
        temperature = self.simulate_thermal_case('test_thermal_heat_flux/test_thermal_heat_flux_3D10N')
        self.assertAlmostEqual(self.etalon_value3, temperature[124])

    def test_thermal_heat_flux_3D8N(self):
        temperature = self.simulate_thermal_case('test_thermal_heat_flux/test_thermal_heat_flux_3D8N')
        self.assertAlmostEqual(self.etalon_value4, temperature[38])

    def test_thermal_heat_flux_3D20N(self):
        temperature = self.simulate_thermal_case('test_thermal_heat_flux/test_thermal_heat_flux_3D20N')
        self.assertAlmostEqual(self.etalon_value4, temperature[213])

    def test_transient_thermal_heat_flux_2D6N(self):
        temperature = self.simulate_thermal_case(
            'test_transient_thermal_heat_flux/test_transient_thermal_heat_flux_2D6N')
        self.assertAlmostEqual(0.46919946397780093, temperature[57])

    def test_transient_thermal_heat_flux_2D10N(self):
        temperature = self.simulate_thermal_case(
            'test_transient_thermal_heat_flux/test_transient_thermal_heat_flux_2D10N')
        self.assertAlmostEqual(0.4674055416030332, temperature[77])

    def test_transient_thermal_heat_flux_2D15N(self):
        temperature = self.simulate_thermal_case(
            'test_transient_thermal_heat_flux/test_transient_thermal_heat_flux_2D15N')
        self.assertAlmostEqual(0.46761403285540487, temperature[97])

    def test_transient_thermal_heat_flux_2D4N(self):
        temperature = self.simulate_thermal_case(
            'test_transient_thermal_heat_flux/test_transient_thermal_heat_flux_2D4N')
        self.assertAlmostEqual(0.12253593527932072, temperature[18])

    def test_transient_thermal_heat_flux_2D8N(self):
        temperature = self.simulate_thermal_case(
            'test_transient_thermal_heat_flux/test_transient_thermal_heat_flux_2D8N')
        self.assertAlmostEqual(0.20716154048406607, temperature[50])

    def test_transient_thermal_heat_flux_2D9N(self):
        temperature = self.simulate_thermal_case(
            'test_transient_thermal_heat_flux/test_transient_thermal_heat_flux_2D9N')
        self.assertAlmostEqual(0.20715104139698065, temperature[63])

    def test_transient_thermal_heat_flux_3D10N(self):
        temperature = self.simulate_thermal_case(
            'test_transient_thermal_heat_flux/test_transient_thermal_heat_flux_3D10N')
        self.assertAlmostEqual(1.2348380064219207, temperature[124])

    def test_thermal_fixed_temperature_2D6N(self):
        temperature = self.simulate_thermal_case('test_thermal_fixed_temperature/test_thermal_fixed_temperature_2D6N')
        self.assertAlmostEqual(13.08528783780587, temperature[57])

    def test_thermal_fixed_temperature_2D10N(self):
        temperature = self.simulate_thermal_case('test_thermal_fixed_temperature/test_thermal_fixed_temperature_2D10N')
        self.assertAlmostEqual(13.08528783780587, temperature[77])

    def test_thermal_fixed_temperature_2D15N(self):
        temperature = self.simulate_thermal_case('test_thermal_fixed_temperature/test_thermal_fixed_temperature_2D15N')
        self.assertAlmostEqual(13.08528783780587, temperature[97])

    def test_thermal_fixed_temperature_2D4N(self):
        temperature = self.simulate_thermal_case('test_thermal_fixed_temperature/test_thermal_fixed_temperature_2D4N')
        self.assertAlmostEqual(4.9497705225, temperature[18])

    def test_thermal_fixed_temperature_2D8N(self):
        temperature = self.simulate_thermal_case('test_thermal_fixed_temperature/test_thermal_fixed_temperature_2D8N')
        self.assertAlmostEqual(4.9497705225, temperature[50])

    def test_thermal_fixed_temperature_2D9N(self):
        temperature = self.simulate_thermal_case('test_thermal_fixed_temperature/test_thermal_fixed_temperature_2D9N')
        self.assertAlmostEqual(4.9497705225, temperature[63])

    def test_thermal_fixed_temperature_3D10N(self):
        temperature = self.simulate_thermal_case('test_thermal_fixed_temperature/test_thermal_fixed_temperature_3D10N')
        self.assertAlmostEqual(16.39151949, temperature[124])

    def test_transient_thermal_fixed_temperature_2D6N(self):
        temperature = self.simulate_thermal_case(
            'test_transient_thermal_fixed_temperature/test_transient_thermal_fixed_temperature_2D6N')
        self.assertAlmostEqual(2.934877024764686, temperature[57])

    def test_transient_thermal_fixed_temperature_2D10N(self):
        temperature = self.simulate_thermal_case(
            'test_transient_thermal_fixed_temperature/test_transient_thermal_fixed_temperature_2D10N')
        self.assertAlmostEqual(2.941952981817704, temperature[77])

    def test_transient_thermal_fixed_temperature_2D15N(self):
        temperature = self.simulate_thermal_case(
            'test_transient_thermal_fixed_temperature/test_transient_thermal_fixed_temperature_2D15N')
        self.assertAlmostEqual(2.942246847876168, temperature[97])

    def test_transient_thermal_fixed_temperature_2D4N(self):
        temperature = self.simulate_thermal_case(
            'test_transient_thermal_fixed_temperature/test_transient_thermal_fixed_temperature_2D4N')
        self.assertAlmostEqual(0.2525670084947606, temperature[18])

    def test_transient_thermal_fixed_temperature_2D8N(self):
        temperature = self.simulate_thermal_case(
            'test_transient_thermal_fixed_temperature/test_transient_thermal_fixed_temperature_2D8N')
        self.assertAlmostEqual(0.37765333581693433, temperature[50])

    def test_transient_thermal_fixed_temperature_2D9N(self):
        temperature = self.simulate_thermal_case(
            'test_transient_thermal_fixed_temperature/test_transient_thermal_fixed_temperature_2D9N')
        self.assertAlmostEqual(0.37766517031304647, temperature[63])

    def test_transient_thermal_fixed_temperature_3D10N(self):
        temperature = self.simulate_thermal_case(
            'test_transient_thermal_fixed_temperature/test_transient_thermal_fixed_temperature_3D10N')
        self.assertAlmostEqual(6.001937486605789, temperature[124])

    def test_thermal_line_element_2D3N(self):
        temperature = self.simulate_thermal_case('test_thermal_line_element/test_thermal_line_element_2D3N')
        self.assertAlmostEqual(self.etalon_value5, temperature[2])

    def test_thermal_line_element_2D4N(self):
        temperature = self.simulate_thermal_case('test_thermal_line_element/test_thermal_line_element_2D4N')
        self.assertAlmostEqual(self.etalon_value5, temperature[2])

    def test_thermal_line_element_2D5N(self):
        temperature = self.simulate_thermal_case('test_thermal_line_element/test_thermal_line_element_2D5N')
        self.assertAlmostEqual(self.etalon_value5, temperature[2])

    def test_thermal_line_element_3D3N(self):
        temperature = self.simulate_thermal_case('test_thermal_line_element/test_thermal_line_element_3D3N')
        self.assertAlmostEqual(self.etalon_value5, temperature[2])

    def test_thermal_point_flux_2D3N(self):
        temperature = self.simulate_thermal_case('test_thermal_heat_flux_line_element/test_thermal_point_flux_2D3N')
        self.assertAlmostEqual(self.etalon_value6, temperature[0])

    def test_thermal_point_flux_2D4N(self):
        temperature = self.simulate_thermal_case('test_thermal_heat_flux_line_element/test_thermal_point_flux_2D4N')
        self.assertAlmostEqual(self.etalon_value6, temperature[0])

    def test_thermal_point_flux_2D5N(self):
        temperature = self.simulate_thermal_case('test_thermal_heat_flux_line_element/test_thermal_point_flux_2D5N')
        self.assertAlmostEqual(self.etalon_value6, temperature[0])

    def test_thermal_point_flux_3D3N(self):
        temperature = self.simulate_thermal_case('test_thermal_heat_flux_line_element/test_thermal_point_flux_3D3N')
        self.assertAlmostEqual(self.etalon_value6, temperature[0])

    def test_thermal_filter_element_2D3N(self):
        temperature = self.simulate_thermal_case('test_thermal_filter_element/test_thermal_filter_element_2D3N')
        self.assertAlmostEqual(34.657035390578656, temperature[79])


if __name__ == '__main__':
    KratosUnittest.main()
