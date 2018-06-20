from KratosMultiphysics import *

class AdjointCPPUnitTests(KratosUnittest.TestCase):

    def setUp(self):
        pass

    def testSensitivityCPP(self):
        Tester.SetVerbosity(Tester.Verbosity.TESTS_OUTPUTS)
        Tester.RunTestSuite("KratosSensitivityTestSuite")

    def tearDown(self):
        pass

if __name__ == '__main__':
    test = AdjointCPPUnitTests()
    test.setUp()
    test.testSensitivityCPP()
    test.tearDown()