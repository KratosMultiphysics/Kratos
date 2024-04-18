import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper

class KratosGeoMechanicsTransientPressureLineElementTests(KratosUnittest.TestCase):
    """
    This class contains benchmark tests which are checked with the regression on a previously obtained value.
    """
    etalon_value1 = -20000.0

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def test_oblique_line_element2D2N(self):
        test_name = 'test_oblique_line_element2D2N'
        file_path = test_helper.get_file_path(os.path.join('test_pressure_line_element', test_name))
        simulation = test_helper.run_kratos(file_path)
        pressure = test_helper.get_water_pressure(simulation)
        self.assertAlmostEqual(self.etalon_value1, pressure[2])

    def test_oblique_line_element2D3N(self):
        test_name = 'test_oblique_line_element2D3N'
        file_path = test_helper.get_file_path(os.path.join('test_pressure_line_element', test_name))
        simulation = test_helper.run_kratos(file_path)
        pressure = test_helper.get_water_pressure(simulation)
        self.assertAlmostEqual(self.etalon_value1, pressure[2])

    def test_oblique_line_element2D4N(self):
        test_name = 'test_oblique_line_element2D4N'
        file_path = test_helper.get_file_path(os.path.join('test_pressure_line_element', test_name))
        simulation = test_helper.run_kratos(file_path)
        pressure = test_helper.get_water_pressure(simulation)
        self.assertAlmostEqual(self.etalon_value1, pressure[2])

    def test_oblique_line_element2D5N(self):
        test_name = 'test_oblique_line_element2D5N'
        file_path = test_helper.get_file_path(os.path.join('test_pressure_line_element', test_name))
        simulation = test_helper.run_kratos(file_path)
        pressure = test_helper.get_water_pressure(simulation)
        self.assertAlmostEqual(self.etalon_value1, pressure[2])

    def test_vertical_line_element2D2N(self):
        test_name = 'test_vertical_line_element2D2N'
        file_path = test_helper.get_file_path(os.path.join('test_pressure_line_element', test_name))
        simulation = test_helper.run_kratos(file_path)
        pressure = test_helper.get_water_pressure(simulation)
        self.assertAlmostEqual(self.etalon_value1, pressure[2])
        
    def test_vertical_line_element2D3N(self):
        test_name = 'test_vertical_line_element2D3N'
        file_path = test_helper.get_file_path(os.path.join('test_pressure_line_element', test_name))
        simulation = test_helper.run_kratos(file_path)
        pressure = test_helper.get_water_pressure(simulation)
        self.assertAlmostEqual(self.etalon_value1, pressure[2])
        
    def test_vertical_line_element2D4N(self):
        test_name = 'test_vertical_line_element2D4N'
        file_path = test_helper.get_file_path(os.path.join('test_pressure_line_element', test_name))
        simulation = test_helper.run_kratos(file_path)
        pressure = test_helper.get_water_pressure(simulation)
        self.assertAlmostEqual(self.etalon_value1, pressure[2])
        
    def test_vertical_line_element2D5N(self):
        test_name = 'test_vertical_line_element2D5N'
        file_path = test_helper.get_file_path(os.path.join('test_pressure_line_element', test_name))
        simulation = test_helper.run_kratos(file_path)
        pressure = test_helper.get_water_pressure(simulation)
        self.assertAlmostEqual(self.etalon_value1, pressure[2])

if __name__ == '__main__':
    KratosUnittest.main()
