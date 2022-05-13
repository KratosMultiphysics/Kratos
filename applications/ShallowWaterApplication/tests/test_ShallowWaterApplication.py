# import Kratos
import KratosMultiphysics as KM

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.KratosUnittest import TestLoader

import argparse

# Small tests
from shallow_water_test_factory import TestShallowWaterElement
from shallow_water_test_factory import TestSemiLagrangianShallowWaterElement
from shallow_water_test_factory import TestConservativeResidualViscosity2D3NElement
from shallow_water_test_factory import TestConservativeGradientJump2D3NElement
from shallow_water_test_factory import TestConservativeFluxCorrected2D3NElement
from shallow_water_test_factory import TestBoussinesq2D3NElement
from shallow_water_test_factory import TestSetTopographyProcess
from shallow_water_test_factory import TestVisualizationMeshProcess
from shallow_water_test_factory import TestMacDonaldShockBenchmark
from shallow_water_test_factory import TestMacDonaldTransitionBenchmark
from shallow_water_test_factory import TestDamBreakBenchmark
from shallow_water_test_factory import TestDryDamBreakBenchmark
from shallow_water_test_factory import TestPlanarSurfaceInParabolaBenchmark
from shallow_water_test_factory import TestSolitaryWaveBenchmark
from shallow_water_test_factory import TestMeshMovingStrategy
from processes_tests.test_line_graph_output_process import TestLineGraphOutputProcess
from processes_tests.test_derivatives_recovery_process import TestDerivativesRecoveryProcess
from processes_tests.test_convergence_output_process import TestConvergenceOutputProcess
from processes_tests.test_wave_generator_process import TestWaveGeneratorProcess

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
    smallSuite.addTests(TestLoader().loadTestsFromTestCases([TestConservativeResidualViscosity2D3NElement]))
    smallSuite.addTests(TestLoader().loadTestsFromTestCases([TestConservativeGradientJump2D3NElement]))
    smallSuite.addTests(TestLoader().loadTestsFromTestCases([TestConservativeFluxCorrected2D3NElement]))
    smallSuite.addTests(TestLoader().loadTestsFromTestCases([TestBoussinesq2D3NElement]))
    smallSuite.addTests(TestLoader().loadTestsFromTestCases([TestSetTopographyProcess]))
    smallSuite.addTests(TestLoader().loadTestsFromTestCases([TestVisualizationMeshProcess]))
    smallSuite.addTests(TestLoader().loadTestsFromTestCases([TestMacDonaldShockBenchmark]))
    smallSuite.addTests(TestLoader().loadTestsFromTestCases([TestMacDonaldTransitionBenchmark]))
    smallSuite.addTests(TestLoader().loadTestsFromTestCases([TestDamBreakBenchmark]))
    smallSuite.addTests(TestLoader().loadTestsFromTestCases([TestDryDamBreakBenchmark]))
    smallSuite.addTests(TestLoader().loadTestsFromTestCases([TestPlanarSurfaceInParabolaBenchmark]))
    smallSuite.addTests(TestLoader().loadTestsFromTestCases([TestSolitaryWaveBenchmark]))
    smallSuite.addTests(TestLoader().loadTestsFromTestCases([TestLineGraphOutputProcess]))
    smallSuite.addTests(TestLoader().loadTestsFromTestCases([TestDerivativesRecoveryProcess]))
    smallSuite.addTests(TestLoader().loadTestsFromTestCases([TestConvergenceOutputProcess]))
    smallSuite.addTests(TestLoader().loadTestsFromTestCases([TestWaveGeneratorProcess]))

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

def run():
    KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
    KratosUnittest.runTests(AssembleTestSuites())

exit_codes = {
    0 : 'OK',
    1 : 'Failed'
}

if __name__ == '__main__':
    parser = argparse.ArgumentParser() # this is a very bad practice, since this overwrites the arguments defined in KratosUnittest
    parser.add_argument('-c', '--cpp', action='store_true')
    args = parser.parse_args()

    KM.Tester.SetVerbosity(KM.Tester.Verbosity.FAILED_TESTS_OUTPUTS)
    cpp_exit_code = KM.Tester.RunTestSuite("ShallowWaterApplicationFastSuite")

    if not args.cpp:
        py_exit_code = run()

    print()
    print('ShallowWaterApplication tests:')
    print('cpp tests ..........', exit_codes[cpp_exit_code])
    if not args.cpp:
        print('python tests .......', exit_codes[py_exit_code])
