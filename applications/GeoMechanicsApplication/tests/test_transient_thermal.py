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
        temp = temperature[49]
        self.assertAlmostEqual(14.17464130639103, temp)

    def test_thermal_fixed_temperature_2D6N(self):
        test_name = 'test_thermal_fixed_temperature_2D6N'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[99]
        self.assertAlmostEqual(6.596927051998985, temp)

    def test_thermal_heat_flux_2D6N(self):
        test_name = 'test_thermal_heat_flux_2D6N'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)
        temperature = test_helper.get_temperature(simulation)
        temp = temperature[224]
        self.assertAlmostEqual(41.93149511165559, temp)


if __name__ == '__main__':
    KratosUnittest.main()
