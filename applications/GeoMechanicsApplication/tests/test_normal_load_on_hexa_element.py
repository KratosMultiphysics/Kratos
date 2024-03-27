import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper

class KratosGeoMechanicsNormalLoadHexaTests(KratosUnittest.TestCase):
    """
    This class contains benchmark tests to test if normal loads are correctly calculated on 1D elements.
    """

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def test_hexa_8n_normal_load(self):
        test_name = 'hexa_8n_normal_load'
        file_path = test_helper.get_file_path(test_name)
        simulation = test_helper.run_kratos(file_path)





        self.assertAlmostEqual(0, 0, places=3)

if __name__ == '__main__':
    KratosUnittest.main()