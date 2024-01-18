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

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def test_thermal_heat_flux_2D3N(self):
        test_name = 'test_thermal_heat_flux_2D3N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', 'test_thermal_heat_flux', test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[37]
        self.assertAlmostEqual(self.etalon_value1, temp)

    def test_thermal_heat_flux_2D6N(self):
        test_name = 'test_thermal_heat_flux_2D6N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', 'test_thermal_heat_flux', test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[57]
        self.assertAlmostEqual(self.etalon_value1, temp)

    def test_thermal_heat_flux_2D10N(self):
        test_name = 'test_thermal_heat_flux_2D10N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', 'test_thermal_heat_flux', test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[77]
        self.assertAlmostEqual(self.etalon_value1, temp)

    def test_thermal_heat_flux_2D15N(self):
        test_name = 'test_thermal_heat_flux_2D15N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', 'test_thermal_heat_flux', test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[97]
        self.assertAlmostEqual(self.etalon_value1, temp)

    def test_thermal_heat_flux_2D4N(self):
        test_name = 'test_thermal_heat_flux_2D4N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', 'test_thermal_heat_flux', test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[18]
        self.assertAlmostEqual(self.etalon_value2, temp)

    def test_thermal_heat_flux_2D8N(self):
        test_name = 'test_thermal_heat_flux_2D8N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', 'test_thermal_heat_flux', test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[50]
        self.assertAlmostEqual(self.etalon_value2, temp)

    def test_thermal_heat_flux_2D9N(self):
        test_name = 'test_thermal_heat_flux_2D9N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', 'test_thermal_heat_flux', test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[63]
        self.assertAlmostEqual(self.etalon_value2, temp)

    def test_thermal_heat_flux_3D4N(self):
        test_name = 'test_thermal_heat_flux_3D4N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', 'test_thermal_heat_flux', test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[22]
        self.assertAlmostEqual(self.etalon_value3, temp)

    def test_thermal_heat_flux_3D10N(self):
        test_name = 'test_thermal_heat_flux_3D10N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', 'test_thermal_heat_flux', test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[124]
        self.assertAlmostEqual(self.etalon_value3, temp)

    def test_transient_thermal_heat_flux_2D3N(self):
        test_name = 'test_transient_thermal_heat_flux_2D3N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', 'test_transient_thermal_heat_flux', test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[37]
        self.assertAlmostEqual(0.3618991346235092, temp)

    def test_transient_thermal_heat_flux_2D6N(self):
        test_name = 'test_transient_thermal_heat_flux_2D6N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', 'test_transient_thermal_heat_flux', test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[57]
        self.assertAlmostEqual(0.46919946397780093, temp)

    def test_transient_thermal_heat_flux_2D10N(self):
        test_name = 'test_transient_thermal_heat_flux_2D10N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', 'test_transient_thermal_heat_flux', test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[77]
        self.assertAlmostEqual(0.4674055416030332, temp)

    def test_transient_thermal_heat_flux_2D15N(self):
        test_name = 'test_transient_thermal_heat_flux_2D15N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', 'test_transient_thermal_heat_flux', test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[97]
        self.assertAlmostEqual(0.46761403285540487, temp)

    def test_transient_thermal_heat_flux_2D4N(self):
        test_name = 'test_transient_thermal_heat_flux_2D4N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', 'test_transient_thermal_heat_flux', test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[18]
        self.assertAlmostEqual(0.12253593527932072, temp)

    def test_transient_thermal_heat_flux_2D8N(self):
        test_name = 'test_transient_thermal_heat_flux_2D8N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', 'test_transient_thermal_heat_flux', test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[50]
        self.assertAlmostEqual(0.20716154048406607, temp)

    def test_transient_thermal_heat_flux_2D9N(self):
        test_name = 'test_transient_thermal_heat_flux_2D9N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', 'test_transient_thermal_heat_flux', test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[63]
        self.assertAlmostEqual(0.20715104139698065, temp)

    def test_transient_thermal_heat_flux_3D4N(self):
        test_name = 'test_transient_thermal_heat_flux_3D4N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', 'test_transient_thermal_heat_flux', test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[22]
        self.assertAlmostEqual(0.8936587648750058, temp)

    def test_transient_thermal_heat_flux_3D10N(self):
        test_name = 'test_transient_thermal_heat_flux_3D10N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', 'test_transient_thermal_heat_flux', test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[124]
        self.assertAlmostEqual(1.2294110493528096, temp)

    def test_thermal_fixed_temperature_2D3N(self):
        test_name = 'test_thermal_fixed_temperature_2D3N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', 'test_thermal_fixed_temperature', test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[37]
        self.assertAlmostEqual(13.08528783780587, temp)

    def test_thermal_fixed_temperature_2D6N(self):
        test_name = 'test_thermal_fixed_temperature_2D6N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', 'test_thermal_fixed_temperature', test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[57]
        self.assertAlmostEqual(13.08528783780587, temp)

    def test_thermal_fixed_temperature_2D10N(self):
        test_name = 'test_thermal_fixed_temperature_2D10N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', 'test_thermal_fixed_temperature', test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[77]
        self.assertAlmostEqual(13.08528783780587, temp)

    def test_thermal_fixed_temperature_2D15N(self):
        test_name = 'test_thermal_fixed_temperature_2D15N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', 'test_thermal_fixed_temperature', test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[97]
        self.assertAlmostEqual(13.08528783780587, temp)

    def test_thermal_fixed_temperature_2D4N(self):
        test_name = 'test_thermal_fixed_temperature_2D4N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', 'test_thermal_fixed_temperature', test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[18]
        self.assertAlmostEqual(4.9497705225, temp)

    def test_thermal_fixed_temperature_2D8N(self):
        test_name = 'test_thermal_fixed_temperature_2D8N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', 'test_thermal_fixed_temperature', test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[50]
        self.assertAlmostEqual(4.9497705225, temp)

    def test_thermal_fixed_temperature_2D9N(self):
        test_name = 'test_thermal_fixed_temperature_2D9N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', 'test_thermal_fixed_temperature', test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[63]
        self.assertAlmostEqual(4.9497705225, temp)

    def test_thermal_fixed_temperature_3D4N(self):
        test_name = 'test_thermal_fixed_temperature_3D4N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', 'test_thermal_fixed_temperature', test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[22]
        self.assertAlmostEqual(16.39151949, temp)

    def test_thermal_fixed_temperature_3D10N(self):
        test_name = 'test_thermal_fixed_temperature_3D10N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', 'test_thermal_fixed_temperature', test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[124]
        self.assertAlmostEqual(16.39151949, temp)

    def test_transient_thermal_fixed_temperature_2D3N(self):
        test_name = 'test_transient_thermal_fixed_temperature_2D3N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', 'test_transient_thermal_fixed_temperature', test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[37]
        self.assertAlmostEqual(3.1312633472490803, temp)

    def test_transient_thermal_fixed_temperature_2D6N(self):
        test_name = 'test_transient_thermal_fixed_temperature_2D6N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', 'test_transient_thermal_fixed_temperature', test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[57]
        self.assertAlmostEqual(2.9280665380753517, temp)

    def test_transient_thermal_fixed_temperature_2D10N(self):
        test_name = 'test_transient_thermal_fixed_temperature_2D10N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', 'test_transient_thermal_fixed_temperature', test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[77]
        self.assertAlmostEqual(2.98104309613254, temp)

    def test_transient_thermal_fixed_temperature_2D15N(self):
        test_name = 'test_transient_thermal_fixed_temperature_2D15N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', 'test_transient_thermal_fixed_temperature', test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[97]
        self.assertAlmostEqual(2.941069197319062, temp)

    def test_transient_thermal_fixed_temperature_2D4N(self):
        test_name = 'test_transient_thermal_fixed_temperature_2D4N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', 'test_transient_thermal_fixed_temperature', test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[18]
        self.assertAlmostEqual(0.36178961457330816, temp)

    def test_transient_thermal_fixed_temperature_2D8N(self):
        test_name = 'test_transient_thermal_fixed_temperature_2D8N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', 'test_transient_thermal_fixed_temperature', test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[50]
        self.assertAlmostEqual(0.3723784708947167, temp)

    def test_transient_thermal_fixed_temperature_2D9N(self):
        test_name = 'test_transient_thermal_fixed_temperature_2D9N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', 'test_transient_thermal_fixed_temperature', test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[63]
        self.assertAlmostEqual(0.37239021552505724, temp)
        
    def test_transient_thermal_fixed_temperature_3D4N(self):
        test_name = 'test_transient_thermal_fixed_temperature_3D4N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', 'test_transient_thermal_fixed_temperature', test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[22]
        self.assertAlmostEqual(7.49001003417586, temp)

    def test_transient_thermal_fixed_temperature_3D10N(self):
        test_name = 'test_transient_thermal_fixed_temperature_3D10N'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', 'test_transient_thermal_fixed_temperature', test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[124]
        self.assertAlmostEqual(5.970566939746188, temp)

    def test_micro_climate_1(self):
        test_name = 'test_micro_climate_1'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[4]
        self.assertAlmostEqual(4.409776948705066, temp)
        
    def test_micro_climate_2(self):
        test_name = 'test_micro_climate_2'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[4]
        self.assertAlmostEqual(4.088147833943762, temp)

    def test_micro_climate_3(self):
        test_name = 'test_micro_climate_3'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[4]
        self.assertAlmostEqual(4.507820382303552, temp)

    def test_micro_climate_4(self):
        test_name = 'test_micro_climate_4'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[4]
        self.assertAlmostEqual(6.213353113038092, temp)

    def test_micro_climate_5(self):
        test_name = 'test_micro_climate_5'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[4]
        self.assertAlmostEqual(4.507820382351035, temp)

    def test_micro_climate_6(self):
        test_name = 'test_micro_climate_6'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[4]
        self.assertAlmostEqual(6.366392882971179, temp)

    def test_micro_climate_7(self):
        test_name = 'test_micro_climate_7'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[4]
        self.assertAlmostEqual(7.1712783573336365, temp)

    def test_micro_climate_8(self):
        test_name = 'test_micro_climate_8'
        file_path = test_helper.get_file_path(os.path.join('test_thermal_element', test_name))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[4]
        self.assertAlmostEqual(6.1263675349643965, temp)

if __name__ == '__main__':
    KratosUnittest.main()
