import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper


class KratosGeoMechanicsTransientThermalTests(KratosUnittest.TestCase):
    """
    This class contains benchmark tests which are checked with the regression on a previously obtained value.
    """
    etalon_value1 = 28.04411163544510063559
    etalon_value2 = 41.3797035928672316
    etalon_value3 = 26.6666666666666667
    etalon_value4 = 88.2768361581920904

    def simulate_thermal_case(self, test_name):
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', test_name))
        simulation = test_helper.run_kratos(file_path)
        return test_helper.get_temperature(simulation)

    def test_thermal_heat_flux_2D3N(self):
        temperature = self.simulate_thermal_case('test_thermal_heat_flux/test_thermal_heat_flux_2D3N')
        self.assertAlmostEqual(self.etalon_value1, temperature[37])

    def test_thermal_heat_flux_3D4N(self):
        temperature = self.simulate_thermal_case('test_thermal_heat_flux/test_thermal_heat_flux_3D4N')
        self.assertAlmostEqual(self.etalon_value2, temperature[22])

    def test_transient_thermal_heat_flux_2D3N(self):
        temperature = self.simulate_thermal_case(
            'test_transient_thermal_heat_flux/test_transient_thermal_heat_flux_2D3N')
        self.assertAlmostEqual(0.3618991346235092, temperature[37])

    def test_transient_thermal_heat_flux_3D4N(self):
        temperature = self.simulate_thermal_case(
            'test_transient_thermal_heat_flux/test_transient_thermal_heat_flux_3D4N')
        self.assertAlmostEqual(0.8936587648750058, temperature[22])

    def test_thermal_fixed_temperature_2D3N(self):
        temperature = self.simulate_thermal_case('test_thermal_fixed_temperature/test_thermal_fixed_temperature_2D3N')
        self.assertAlmostEqual(13.08528783780587, temperature[37])

    def test_thermal_fixed_temperature_2D3N_newmark(self):
        temperature = self.simulate_thermal_case(
            'test_thermal_fixed_temperature/test_thermal_fixed_temperature_2D3N_newmark')
        self.assertAlmostEqual(13.08528783780587, temperature[37])

    def test_thermal_fixed_temperature_3D4N(self):
        temperature = self.simulate_thermal_case('test_thermal_fixed_temperature/test_thermal_fixed_temperature_3D4N')
        self.assertAlmostEqual(16.39151949, temperature[22])

    def test_transient_thermal_fixed_temperature_2D3N(self):
        temperature = self.simulate_thermal_case(
            'test_transient_thermal_fixed_temperature/test_transient_thermal_fixed_temperature_2D3N')
        self.assertAlmostEqual(2.711294346531134, temperature[37])

    def test_transient_thermal_fixed_temperature_3D4N(self):
        temperature = self.simulate_thermal_case(
            'test_transient_thermal_fixed_temperature/test_transient_thermal_fixed_temperature_3D4N')
        self.assertAlmostEqual(6.368655984441822, temperature[22])

    def test_micro_climate_1(self):
        temperature = self.simulate_thermal_case('test_micro_climate_1')
        self.assertAlmostEqual(4.409776948705066, temperature[4])

    def test_micro_climate_2(self):
        temperature = self.simulate_thermal_case('test_micro_climate_2')
        self.assertAlmostEqual(4.088147833943762, temperature[4])

    def test_micro_climate_3(self):
        temperature = self.simulate_thermal_case('test_micro_climate_3')
        self.assertAlmostEqual(4.507820382303552, temperature[4])

    def test_micro_climate_4(self):
        temperature = self.simulate_thermal_case('test_micro_climate_4')
        self.assertAlmostEqual(6.213353113038092, temperature[4])

    def test_micro_climate_5(self):
        temperature = self.simulate_thermal_case('test_micro_climate_5')
        self.assertAlmostEqual(4.507820382351035, temperature[4])

    def test_micro_climate_6(self):
        temperature = self.simulate_thermal_case('test_micro_climate_6')
        self.assertAlmostEqual(6.366392882971179, temperature[4])

    def test_micro_climate_7(self):
        temperature = self.simulate_thermal_case('test_micro_climate_7')
        self.assertAlmostEqual(7.1712783573336365, temperature[4])

    def test_micro_climate_8(self):
        temperature = self.simulate_thermal_case('test_micro_climate_8')
        self.assertAlmostEqual(6.1263675349643965, temperature[4])

    def test_thermal_line_element_2D2N(self):
        temperature = self.simulate_thermal_case('test_thermal_line_element/test_thermal_line_element_2D2N')
        self.assertAlmostEqual(self.etalon_value3, temperature[2])

    def test_thermal_line_element_3D2N(self):
        temperature = self.simulate_thermal_case('test_thermal_line_element/test_thermal_line_element_3D2N')
        self.assertAlmostEqual(self.etalon_value3, temperature[2])

    def test_thermal_point_flux_2D2N(self):
        temperature = self.simulate_thermal_case('test_thermal_heat_flux_line_element/test_thermal_point_flux_2D2N')
        self.assertAlmostEqual(self.etalon_value4, temperature[0])

    def test_thermal_point_flux_3D2N(self):
        temperature = self.simulate_thermal_case('test_thermal_heat_flux_line_element/test_thermal_point_flux_3D2N')
        self.assertAlmostEqual(self.etalon_value4, temperature[0])

    def test_thermal_filter_element_2D2N(self):
        temperature = self.simulate_thermal_case('test_thermal_filter_element/test_thermal_filter_element_2D2N')
        self.assertAlmostEqual(34.65690605787046, temperature[22])


if __name__ == '__main__':
    KratosUnittest.main()
