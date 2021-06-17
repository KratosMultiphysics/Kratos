# import Kratos
import KratosMultiphysics as KM

# cpp tests
import run_cpp_unit_tests

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Small tests
from shallow_water_test_factory import TestShallowWaterElement
from shallow_water_test_factory import TestSemiLagrangianShallowWaterElement
from shallow_water_test_factory import TestShallowWater2D3NElement
from shallow_water_test_factory import TestMonotonicShallowWater2D3NElement
from shallow_water_test_factory import TestSetTopographyProcess
from shallow_water_test_factory import TestVisualizationMeshProcess
from shallow_water_test_factory import TestNodesOutputProcess
from shallow_water_test_factory import TestMacDonaldShockBenchmark
from shallow_water_test_factory import TestDamBreakBenchmark
from shallow_water_test_factory import TestDryDamBreakBenchmark
from shallow_water_test_factory import TestPlanarSurfaceInParabolaBenchmark
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
    smallSuite.addTest(TestShallowWater2D3NElement('test_execution'))
    smallSuite.addTest(TestMonotonicShallowWater2D3NElement('test_execution'))
    smallSuite.addTest(TestSetTopographyProcess('test_execution'))
    smallSuite.addTest(TestVisualizationMeshProcess('test_execution'))
    smallSuite.addTest(TestNodesOutputProcess('test_execution'))
    smallSuite.addTest(TestMacDonaldShockBenchmark('test_execution'))
    smallSuite.addTest(TestDamBreakBenchmark('test_execution'))
    smallSuite.addTest(TestDryDamBreakBenchmark('test_execution'))
    smallSuite.addTest(TestPlanarSurfaceInParabolaBenchmark('test_execution'))
    smallSuite.addTests(_loadTestsFromTestCases(TestConvergenceOutputProcess))

    # Create a test suit with the selected tests plus all small tests
    nightlySuite = suites['nightly']
    nightlySuite.addTests(smallSuite)
    nightlySuite.addTest(TestShallowWaterElement('test_execution'))
    nightlySuite.addTest(TestSemiLagrangianShallowWaterElement('test_execution'))

    # Create a test suit that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(nightlySuite)

    return suites

def _loadTestsFromTestCases(test_case_module):
    return KratosUnittest.TestLoader().loadTestsFromTestCases([test_case_module])

if __name__ == '__main__':
    KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
    run_cpp_unit_tests.run()
    KratosUnittest.runTests(AssembleTestSuites())
