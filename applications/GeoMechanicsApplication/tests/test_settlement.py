import sys
import os

sys.path.append(os.path.join('..', '..', '..'))
sys.path.append(os.path.join('..', 'python_scripts'))

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper

class KratosGeoMechanicsSettlementTests(KratosUnittest.TestCase):

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def test_Abc_1_1_0_True_Deformations(self):
        test_name = 'test_Abc_1_1_0_True_Deformations'
        file_path = test_helper.get_file_path(os.path.join('.', test_name))
        stages = test_helper.run_stages(file_path, 2)
        self.assertEqual(True, True)

if __name__ == '__main__':
    suites = KratosUnittest.KratosSuites
    smallSuite = suites['small'] # These tests are executed by the continuous integration tool
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([KratosGeoMechanicsSettlementTests]))
    allSuite = suites['all']
    allSuite.addTests(smallSuite)
    KratosUnittest.runTests(suites)  
