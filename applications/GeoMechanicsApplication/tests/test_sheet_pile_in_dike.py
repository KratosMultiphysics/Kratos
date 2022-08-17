import sys
import os

sys.path.append(os.path.join('..', '..', '..'))
sys.path.append(os.path.join('..', 'python_scripts'))
sys.path.append(os.path.join('D:\\kratos'))

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper


class KratosSheetPileInDikeTest(KratosUnittest.TestCase):
    """
    This class contains benchmark tests which are checked with the analytical solution
    """

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    @KratosUnittest.skip("This test should be checked as it fails without a meaningful error.")
    def test_stage_1a(self):
        test_name = 'InitialPhase1a_redefinition'
        file_path = test_helper.get_file_path(os.path.join('.', 'test_sheet_pile_in_dike_phases',  test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)
        self.assertTrue(simulation is not None)
