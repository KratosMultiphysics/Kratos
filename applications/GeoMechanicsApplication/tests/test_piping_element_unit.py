import sys
import os

# sys.path.append(os.path.join('..', '..', '..'))
# sys.path.append(os.path.join('..', 'python_scripts'))
# sys.path.append(r"D:\software_development\Kratos\bin\Debug")

from KratosMultiphysics import Tester
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.GeoMechanicsApplication

class TestUnitPipingElements(KratosUnittest.TestCase):

    # Tester.RunAllTestCases() #Test all cases
    def setUp(self):
        Tester.SetVerbosity(Tester.Verbosity.TESTS_OUTPUTS) # Set the verbosity level
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def test_water_pressure_gradient(self):
        exitcode = Tester.RunTestCases("TestCalculateWaterPressureGradient")
        self.assertTrue(exitcode == 0)

    def test_equilibrium_pipe_height(self):
        exitcode = Tester.RunTestCases("TestCalculateEquilibriumPipeHeight")
        self.assertTrue(exitcode == 0)



if __name__ == '__main__':
    # Tester.RunTestSuite("KratosGeoMechanicsFastSuite")
    # Tester.RunAllTestCases()
    KratosUnittest.main()
