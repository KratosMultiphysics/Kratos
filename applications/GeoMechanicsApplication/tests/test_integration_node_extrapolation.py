import os
from functools import reduce
import math

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper





class TestExtrapolation(KratosUnittest.TestCase):
    """
    Class that contains elementary groundwater flow tests.

    """

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def test_3_noded_triangle(self):
        """ Triangle2D3N """
        test_name = 'test_integration_node_extrapolation//3_noded_triangle'
        file_path = test_helper.get_file_path(os.path.join('.', test_name))
        simulation = test_helper.run_kratos(file_path)
       
    def test_6_noded_triangle(self):
        """ Triangle2D6N """
        test_name = 'test_integration_node_extrapolation//6_noded_triangle'
        file_path = test_helper.get_file_path(os.path.join('.', test_name))
        simulation = test_helper.run_kratos(file_path)
        
if __name__=="__main__":
    tests = TestExtrapolation()
    tests.test_3_noded_triangle()
    tests.test_6_noded_triangle()