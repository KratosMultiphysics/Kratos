from KratosMultiphysics import *
from KratosMultiphysics.DelaunayMeshingApplication import *

import KratosMultiphysics.KratosUnittest as KratosUnittest


def GetFilePath(fileName):
    return os.path.dirname(__file__) + "/" + fileName


class KratosDelaunayMeshingGeneralTests(KratosUnittest.TestCase):

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def testSmallExample(self):
        self.AssertEqual(True, True)

    def testNightlyFirstExample(self):
        self.AssertEqual(True, True)

    def testNightlySecondExample(self):
        self.AssertEqual(True, True)
