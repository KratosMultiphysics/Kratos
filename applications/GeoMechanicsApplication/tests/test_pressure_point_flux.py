import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper

class KratosGeoMechanicsTransientPressurePointFluxTests(KratosUnittest.TestCase):
    """
    This class contains benchmark tests which are checked with the regression on a previously obtained value.
    """
    etalon_value1 = -21000.0

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def check_water_pressure(self, test_name, etalon_value):
        file_path = test_helper.get_file_path(os.path.join('test_pressure_point_flux', test_name))
        simulation = test_helper.run_kratos(file_path)
        pressure = test_helper.get_water_pressure(simulation)
        self.assertAlmostEqual(etalon_value, pressure[2])

    def test_pressure_point_flux_line_element2D2N(self):
        self.check_water_pressure("test_pressure_point_flux_line_element2D2N", self.etalon_value1)
        
    def test_pressure_point_flux_line_element3D2N(self):
        self.check_water_pressure("test_pressure_point_flux_line_element3D2N", self.etalon_value1)

if __name__ == '__main__':
    KratosUnittest.main()
