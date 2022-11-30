import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper


class TestKRATOSGEO_14_hydrostatic_case(KratosUnittest.TestCase):

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def test_1(self):
        test_name = 'test_bugfixes/KRATOSGEO_14_hydrostatic_case/test_1.gid'
        file_path = test_helper.get_file_path(os.path.join('.', test_name))
        simulation = test_helper.run_kratos(file_path)
        self.assertTrue(abs(min(test_helper.get_water_pressure(simulation)) - (-20) ) < 1e-5)
        self.assertTrue(abs(max(test_helper.get_water_pressure(simulation)) - 0) < 1e-5)

    def test_2(self):
        test_name = 'test_bugfixes/KRATOSGEO_14_hydrostatic_case/test_2.gid'
        file_path = test_helper.get_file_path(os.path.join('.', test_name))
        simulation = test_helper.run_kratos(file_path)
        self.assertTrue(abs(min(test_helper.get_water_pressure(simulation)) - (-20) ) < 1e-5)
        self.assertTrue(abs(max(test_helper.get_water_pressure(simulation)) - 0) < 1e-5)

    def test_3(self):
        test_name = 'test_bugfixes/KRATOSGEO_14_hydrostatic_case/test_3.gid'
        file_path = test_helper.get_file_path(os.path.join('.', test_name))
        simulation = test_helper.run_kratos(file_path)
        self.assertTrue(abs(min(test_helper.get_water_pressure(simulation)) - (-20) ) < 1e-5)
        self.assertTrue(abs(max(test_helper.get_water_pressure(simulation)) - 0) < 1e-5)

    def test_4(self):
        test_name = 'test_bugfixes/KRATOSGEO_14_hydrostatic_case/test_4.gid'
        file_path = test_helper.get_file_path(os.path.join('.', test_name))
        simulation = test_helper.run_kratos(file_path)
        self.assertTrue(abs(min(test_helper.get_water_pressure(simulation)) - (-20) ) < 1e-5)
        self.assertTrue(abs(max(test_helper.get_water_pressure(simulation)) - 0) < 1e-5)


if __name__ == '__main__':
    KratosUnittest.main()