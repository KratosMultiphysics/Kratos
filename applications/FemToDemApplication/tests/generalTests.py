import KratosMultiphysics
import KratosMultiphysics.FemToDemApplication

import KratosMultiphysics.KratosUnittest as KratosUnittest


def GetFilePath(fileName):
    return os.path.dirname(__file__) + "/" + fileName


class KratosFemToDemGeneralTests(KratosUnittest.TestCase):

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def testSmallExample(self):
        self.assertEqual(True, True)

    def testNightlyFirstExample(self):
        self.assertEqual(True, True)

    def testNightlySecondExample(self):
        self.assertEqual(True, True)
