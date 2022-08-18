from KratosMultiphysics import Tester
import os
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestUnitPipingElements(KratosUnittest.TestCase):

    def setUp(self):
        Tester.SetVerbosity(Tester.Verbosity.TESTS_OUTPUTS) # Set the verbosity level

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def test_water_pressure_gradient(self):
        print(os.getcwd())
        exitcode = Tester.RunTestCases("TestCalculateWaterPressureGradient")
        self.assertTrue(exitcode == 0)

    def test_equilibrium_pipe_height(self):
        print(os.getcwd())
        exitcode = Tester.RunTestCases("TestCalculateEquilibriumPipeHeight")
        self.assertTrue(exitcode == 0)

    def test_ErosionProcessStrategy(self):
        print(os.getcwd())
        exitcode = Tester.RunTestCases("TestErosionProcessStrategy")
        self.assertTrue(exitcode == 0)


if __name__ == '__main__':
    KratosUnittest.main()
