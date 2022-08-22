from KratosMultiphysics import Tester
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestNormalFluxCondition(KratosUnittest.TestCase):

    def setUp(self):
        Tester.SetVerbosity(Tester.Verbosity.TESTS_OUTPUTS) # Set the verbosity level

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def test_horizontal_normal_flux(self):
        exitcode = Tester.RunTestCases("TestCalculateHorizontalNormalFlux")
        self.assertTrue(exitcode == 0)

    def test_inclined_normal_flux(self):
        exitcode = Tester.RunTestCases("TestCalculateInclinedNormalFlux")
        self.assertTrue(exitcode == 0)

if __name__ == '__main__':
    KratosUnittest.main()