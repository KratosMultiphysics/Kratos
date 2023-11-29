import sys
import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper

class KratosGeoMechanicsTransientThermalTests(KratosUnittest.TestCase):
    """
    This class contains benchmark tests which are checked with the analytical solution
    """

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def test_thermal_fixed_temperature_2D3N(self):
        test_name = 'test_thermal_fixed_temperature_2D3N'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[37]
        self.assertAlmostEqual(13.08528783780587, temp)

    def test_thermal_fixed_temperature_2D6N(self):
        test_name = 'test_thermal_fixed_temperature_2D6N'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[57]
        self.assertAlmostEqual(13.08528783780587, temp)

    def test_thermal_fixed_temperature_2D10N(self):
        test_name = 'test_thermal_fixed_temperature_2D10N'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[77]
        self.assertAlmostEqual(13.08528783780587, temp)

    def test_thermal_fixed_temperature_2D15N(self):
        test_name = 'test_thermal_fixed_temperature_2D15N'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[97]
        self.assertAlmostEqual(13.08528783780587, temp)

    def test_thermal_fixed_temperature_2D4N(self):
        test_name = 'test_thermal_fixed_temperature_2D4N'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[18]
        self.assertAlmostEqual(4.9497705225, temp)
        
    def test_thermal_fixed_temperature_2D8N(self):
        test_name = 'test_thermal_fixed_temperature_2D8N'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[50]
        self.assertAlmostEqual(4.9497705225, temp)
        
    def test_thermal_fixed_temperature_2D9N(self):
        test_name = 'test_thermal_fixed_temperature_2D9N'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[63]
        self.assertAlmostEqual(4.9497705225, temp)

    def test_thermal_heat_flux_2D3N(self):
        test_name = 'test_thermal_heat_flux_2D3N'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[37]
        self.assertAlmostEqual(28.04411163544510063559, temp)
        
    def test_thermal_heat_flux_2D6N(self):
        test_name = 'test_thermal_heat_flux_2D6N'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[57]
        self.assertAlmostEqual(28.04411163544510063559, temp)
        
    def test_thermal_heat_flux_2D10N(self):
        test_name = 'test_thermal_heat_flux_2D10N'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[77]
        self.assertAlmostEqual(28.04411163544510063559, temp)
        
    def test_thermal_heat_flux_2D15N(self):
        test_name = 'test_thermal_heat_flux_2D15N'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[97]
        self.assertAlmostEqual(28.04411163544510063559, temp)

    def test_thermal_heat_flux_2D4N(self):
        test_name = 'test_thermal_heat_flux_2D4N'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[18]
        self.assertAlmostEqual(17.55892791313559322, temp)
        
    def test_thermal_heat_flux_2D8N(self):
        test_name = 'test_thermal_heat_flux_2D8N'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[50]
        self.assertAlmostEqual(17.55892791313559322, temp)
        
    def test_thermal_heat_flux_2D9N(self):
        test_name = 'test_thermal_heat_flux_2D9N'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[63]
        self.assertAlmostEqual(17.55892791313559322, temp)
        
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
