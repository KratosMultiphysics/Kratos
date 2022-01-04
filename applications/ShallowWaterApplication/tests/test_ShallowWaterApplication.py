# import Kratos
import KratosMultiphysics as KM

# cpp tests
import run_cpp_unit_tests

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.KratosUnittest import TestLoader

# Small tests
from shallow_water_test_factory import TestShallowWaterElement
from shallow_water_test_factory import TestSemiLagrangianShallowWaterElement
from shallow_water_test_factory import TestShallowWater2D3NElement
from shallow_water_test_factory import TestMonotonicShallowWater2D3NElement
from shallow_water_test_factory import TestSetTopographyProcess
from shallow_water_test_factory import TestVisualizationMeshProcess
from shallow_water_test_factory import TestNodesOutputProcess
from shallow_water_test_factory import TestMacDonaldShockBenchmark
from shallow_water_test_factory import TestMacDonaldTransitionBenchmark
from shallow_water_test_factory import TestDamBreakBenchmark
from shallow_water_test_factory import TestDryDamBreakBenchmark
from shallow_water_test_factory import TestPlanarSurfaceInParabolaBenchmark
from shallow_water_test_factory import TestMeshMovingStrategy
from processes_tests.test_convergence_output_process import TestConvergenceOutputProcess

def AssembleTestSuites():
    ''' Populates the test suites to run.

    Populates the test suites to run. At least, it should populate the suites:
    "small", "nightly" and "all"

    Return
    ------

    suites: A dictionary of suites
        The set of suites with its test_cases added.
    '''
    suites = KratosUnittest.KratosSuites

    # Create a test suit with the selected tests (Small tests):
    smallSuite = suites['small']
    smallSuite.addTests(TestLoader().loadTestsFromTestCases([TestShallowWater2D3NElement]))
    smallSuite.addTests(TestLoader().loadTestsFromTestCases([TestMonotonicShallowWater2D3NElement]))
    smallSuite.addTests(TestLoader().loadTestsFromTestCases([TestSetTopographyProcess]))
    smallSuite.addTests(TestLoader().loadTestsFromTestCases([TestVisualizationMeshProcess]))
    smallSuite.addTests(TestLoader().loadTestsFromTestCases([TestNodesOutputProcess]))
    smallSuite.addTests(TestLoader().loadTestsFromTestCases([TestMacDonaldShockBenchmark]))
    smallSuite.addTests(TestLoader().loadTestsFromTestCases([TestMacDonaldTransitionBenchmark]))
    smallSuite.addTests(TestLoader().loadTestsFromTestCases([TestDamBreakBenchmark]))
    smallSuite.addTests(TestLoader().loadTestsFromTestCases([TestDryDamBreakBenchmark]))
    smallSuite.addTests(TestLoader().loadTestsFromTestCases([TestPlanarSurfaceInParabolaBenchmark]))
    smallSuite.addTests(TestLoader().loadTestsFromTestCases([TestConvergenceOutputProcess]))

    # Create a test suit with the selected tests plus all small tests
    nightlySuite = suites['nightly']
    nightlySuite.addTests(smallSuite)
    nightlySuite.addTests(TestLoader().loadTestsFromTestCases([TestShallowWaterElement]))
    nightlySuite.addTests(TestLoader().loadTestsFromTestCases([TestSemiLagrangianShallowWaterElement]))
    nightlySuite.addTests(TestLoader().loadTestsFromTestCases([TestMeshMovingStrategy]))

    # Create a test suit that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(nightlySuite)

    return suites

if __name__ == '__main__':
    KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
    run_cpp_unit_tests.run()
    KratosUnittest.runTests(AssembleTestSuites())
