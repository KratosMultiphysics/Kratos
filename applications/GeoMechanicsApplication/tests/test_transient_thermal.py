import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper

class KratosGeoMechanicsTransientThermalTests(KratosUnittest.TestCase):
    """
    This class contains benchmark tests which are checked with the regression on a previously obtained value.
    """
    etalon_value1 = 28.04411163544510063559
    etalon_value2 = 17.55892791313559322
    etalon_value3 = 41.3797035928672316
    etalon_value4 = 35.31073446327683615819
    etalon_value5 = 26.6666666666666667
    etalon_value6 = 88.2768361581920904

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass
        
    def check_temperature(self, test_directory, test_name, test_node, etalon_value):
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', test_directory, test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        self.assertAlmostEqual(etalon_value, temperature[test_node])

    def test_thermal_heat_flux_2D3N(self):
        self.check_temperature('test_thermal_heat_flux', 'test_thermal_heat_flux_2D3N', 37, self.etalon_value1)

    def test_thermal_heat_flux_2D6N(self):
        self.check_temperature('test_thermal_heat_flux', 'test_thermal_heat_flux_2D6N', 57, self.etalon_value1)

    def test_thermal_heat_flux_2D10N(self):
        self.check_temperature('test_thermal_heat_flux', 'test_thermal_heat_flux_2D10N', 77, self.etalon_value1)

    def test_thermal_heat_flux_2D15N(self):
        self.check_temperature('test_thermal_heat_flux', 'test_thermal_heat_flux_2D15N', 97, self.etalon_value1)

    def test_thermal_heat_flux_2D4N(self):
        self.check_temperature('test_thermal_heat_flux', 'test_thermal_heat_flux_2D4N', 18, self.etalon_value2)

    def test_thermal_heat_flux_2D8N(self):
        self.check_temperature('test_thermal_heat_flux', 'test_thermal_heat_flux_2D8N', 50, self.etalon_value2)

    def test_thermal_heat_flux_2D9N(self):
        self.check_temperature('test_thermal_heat_flux', 'test_thermal_heat_flux_2D9N', 63, self.etalon_value2)

    def test_thermal_heat_flux_3D4N(self):
        self.check_temperature('test_thermal_heat_flux', 'test_thermal_heat_flux_3D4N', 22, self.etalon_value3)

    def test_thermal_heat_flux_3D10N(self):
        self.check_temperature('test_thermal_heat_flux', 'test_thermal_heat_flux_3D10N', 124, self.etalon_value3)

    def test_thermal_heat_flux_3D8N(self):
        self.check_temperature('test_thermal_heat_flux', 'test_thermal_heat_flux_3D8N', 38, self.etalon_value4)

    def test_thermal_heat_flux_3D20N(self):
        self.check_temperature('test_thermal_heat_flux', 'test_thermal_heat_flux_3D20N', 213, self.etalon_value4)

    def test_transient_thermal_heat_flux_2D3N(self):
        self.check_temperature('test_transient_thermal_heat_flux', 'test_transient_thermal_heat_flux_2D3N', 37, 0.3618991346235092)

    def test_transient_thermal_heat_flux_2D6N(self):
        self.check_temperature('test_transient_thermal_heat_flux', 'test_transient_thermal_heat_flux_2D6N', 57, 0.46919946397780093)

    def test_transient_thermal_heat_flux_2D10N(self):
        self.check_temperature('test_transient_thermal_heat_flux', 'test_transient_thermal_heat_flux_2D10N', 77, 0.4674055416030332)

    def test_transient_thermal_heat_flux_2D15N(self):
        self.check_temperature('test_transient_thermal_heat_flux', 'test_transient_thermal_heat_flux_2D15N', 97, 0.46761403285540487)

    def test_transient_thermal_heat_flux_2D4N(self):
        self.check_temperature('test_transient_thermal_heat_flux', 'test_transient_thermal_heat_flux_2D4N', 18, 0.12253593527932072)

    def test_transient_thermal_heat_flux_2D8N(self):
        self.check_temperature('test_transient_thermal_heat_flux', 'test_transient_thermal_heat_flux_2D8N', 50, 0.20716154048406607)

    def test_transient_thermal_heat_flux_2D9N(self):
        self.check_temperature('test_transient_thermal_heat_flux', 'test_transient_thermal_heat_flux_2D9N', 63, 0.20715104139698065)

    def test_transient_thermal_heat_flux_3D4N(self):
        self.check_temperature('test_transient_thermal_heat_flux', 'test_transient_thermal_heat_flux_3D4N', 22, 0.8936587648750058)

    def test_transient_thermal_heat_flux_3D10N(self):
        self.check_temperature('test_transient_thermal_heat_flux', 'test_transient_thermal_heat_flux_3D10N', 124, 1.2294110493528096)

    def test_thermal_fixed_temperature_2D3N(self):
        self.check_temperature('test_thermal_fixed_temperature', 'test_thermal_fixed_temperature_2D3N', 37, 13.08528783780587)

    def test_thermal_fixed_temperature_2D3N_newmark(self):
        self.check_temperature('test_thermal_fixed_temperature', 'test_thermal_fixed_temperature_2D3N_newmark', 37, 13.08528783780587)

    def test_thermal_fixed_temperature_2D6N(self):
        self.check_temperature('test_thermal_fixed_temperature', 'test_thermal_fixed_temperature_2D6N', 57, 13.08528783780587)

    def test_thermal_fixed_temperature_2D10N(self):
        self.check_temperature('test_thermal_fixed_temperature', 'test_thermal_fixed_temperature_2D10N', 77, 13.08528783780587)

    def test_thermal_fixed_temperature_2D15N(self):
        self.check_temperature('test_thermal_fixed_temperature', 'test_thermal_fixed_temperature_2D15N', 97, 13.08528783780587)

    def test_thermal_fixed_temperature_2D4N(self):
        self.check_temperature('test_thermal_fixed_temperature', 'test_thermal_fixed_temperature_2D4N', 18, 4.9497705225)

    def test_thermal_fixed_temperature_2D8N(self):
        self.check_temperature('test_thermal_fixed_temperature', 'test_thermal_fixed_temperature_2D8N', 50, 4.9497705225)

    def test_thermal_fixed_temperature_2D9N(self):
        self.check_temperature('test_thermal_fixed_temperature', 'test_thermal_fixed_temperature_2D9N', 63, 4.9497705225)

    def test_thermal_fixed_temperature_3D4N(self):
        self.check_temperature('test_thermal_fixed_temperature', 'test_thermal_fixed_temperature_3D4N', 22, 16.39151949)

    def test_thermal_fixed_temperature_3D10N(self):
        self.check_temperature('test_thermal_fixed_temperature', 'test_thermal_fixed_temperature_3D10N', 124, 16.39151949)

    def test_transient_thermal_fixed_temperature_2D3N(self):
        self.check_temperature('test_transient_thermal_fixed_temperature', 'test_transient_thermal_fixed_temperature_2D3N', 37, 3.1312633472490803)

    def test_transient_thermal_fixed_temperature_2D6N(self):
        self.check_temperature('test_transient_thermal_fixed_temperature', 'test_transient_thermal_fixed_temperature_2D6N', 57, 2.9280665380753517)

    def test_transient_thermal_fixed_temperature_2D10N(self):
        self.check_temperature('test_transient_thermal_fixed_temperature', 'test_transient_thermal_fixed_temperature_2D10N', 77, 2.98104309613254)

    def test_transient_thermal_fixed_temperature_2D15N(self):
        self.check_temperature('test_transient_thermal_fixed_temperature', 'test_transient_thermal_fixed_temperature_2D15N', 97, 2.941069197319062)

    def test_transient_thermal_fixed_temperature_2D4N(self):
        self.check_temperature('test_transient_thermal_fixed_temperature', 'test_transient_thermal_fixed_temperature_2D4N', 18, 0.36178961457330816)

    def test_transient_thermal_fixed_temperature_2D8N(self):
        self.check_temperature('test_transient_thermal_fixed_temperature', 'test_transient_thermal_fixed_temperature_2D8N', 50, 0.3723784708947167)

    def test_transient_thermal_fixed_temperature_2D9N(self):
        self.check_temperature('test_transient_thermal_fixed_temperature', 'test_transient_thermal_fixed_temperature_2D9N', 63, 0.37239021552505724)

    def test_transient_thermal_fixed_temperature_3D4N(self):
        self.check_temperature('test_transient_thermal_fixed_temperature', 'test_transient_thermal_fixed_temperature_3D4N', 22, 7.49001003417586)

    def test_transient_thermal_fixed_temperature_3D10N(self):
        self.check_temperature('test_transient_thermal_fixed_temperature', 'test_transient_thermal_fixed_temperature_3D10N', 124, 5.970566939746188)

    def test_micro_climate_1(self):
        self.check_temperature('', 'test_micro_climate_1', 4, 4.409776948705066)

    def test_micro_climate_2(self):
        self.check_temperature('', 'test_micro_climate_2', 4, 4.088147833943762)

    def test_micro_climate_3(self):
        self.check_temperature('', 'test_micro_climate_3', 4, 4.507820382303552)

    def test_micro_climate_4(self):
        self.check_temperature('', 'test_micro_climate_4', 4, 6.213353113038092)

    def test_micro_climate_5(self):
        self.check_temperature('', 'test_micro_climate_5', 4, 4.507820382351035)

    def test_micro_climate_6(self):
        self.check_temperature('', 'test_micro_climate_6', 4, 6.366392882971179)

    def test_micro_climate_7(self):
        self.check_temperature('', 'test_micro_climate_7', 4, 7.1712783573336365)

    def test_micro_climate_8(self):
        self.check_temperature('', 'test_micro_climate_8', 4, 6.1263675349643965)
        
    def test_thermal_line_element_2D2N(self):
        self.check_temperature('test_thermal_line_element', 'test_thermal_line_element_2D2N', 2, self.etalon_value5)

    def test_thermal_line_element_2D3N(self):
        self.check_temperature('test_thermal_line_element', 'test_thermal_line_element_2D3N', 2, self.etalon_value5)

    def test_thermal_line_element_2D4N(self):
        self.check_temperature('test_thermal_line_element', 'test_thermal_line_element_2D4N', 2, self.etalon_value5)

    def test_thermal_line_element_2D5N(self):
        self.check_temperature('test_thermal_line_element', 'test_thermal_line_element_2D5N', 2, self.etalon_value5)

    def test_thermal_line_element_3D2N(self):
        self.check_temperature('test_thermal_line_element', 'test_thermal_line_element_3D2N', 2, self.etalon_value5)

    def test_thermal_line_element_3D3N(self):
        self.check_temperature('test_thermal_line_element', 'test_thermal_line_element_3D3N', 2, self.etalon_value5)

    def test_thermal_point_flux_2D2N(self):
        self.check_temperature('test_thermal_heat_flux_line_element', 'test_thermal_point_flux_2D2N', 0, self.etalon_value6)

    def test_thermal_point_flux_2D3N(self):
        self.check_temperature('test_thermal_heat_flux_line_element', 'test_thermal_point_flux_2D3N', 0, self.etalon_value6)

    def test_thermal_point_flux_2D4N(self):
        self.check_temperature('test_thermal_heat_flux_line_element', 'test_thermal_point_flux_2D4N', 0, self.etalon_value6)

    def test_thermal_point_flux_2D5N(self):
        self.check_temperature('test_thermal_heat_flux_line_element', 'test_thermal_point_flux_2D5N', 0, self.etalon_value6)

    def test_thermal_point_flux_3D2N(self):
        self.check_temperature('test_thermal_heat_flux_line_element', 'test_thermal_point_flux_3D2N', 0, self.etalon_value6)

    def test_thermal_point_flux_3D3N(self):
        self.check_temperature('test_thermal_heat_flux_line_element', 'test_thermal_point_flux_3D3N', 0, self.etalon_value6)

if __name__ == '__main__':
    KratosUnittest.main()
