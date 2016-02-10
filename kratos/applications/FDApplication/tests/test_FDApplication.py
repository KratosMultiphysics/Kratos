import unittest

from KratosMultiphysics import *
from KratosMultiphysics.FDApplication import *


class TestFDApplication(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testNormal(self):
        pass


if __name__ == '__main__':

    # Would run everything
    # unittest.main()

    # would run only test from suit N
    FDApplicationSuit = unittest.TestLoader().loadTestsFromTestCase(TestFDApplication)

    unittest.TextTestRunner(verbosity=2).run(FDApplicationSuit)
