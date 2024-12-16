# import Kratos
import KratosMultiphysics as KM

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.KratosUnittest import TestLoader

# Small tests
from shallow_water_test_factory import TestConservativeResidualViscosity2D3NElement
from shallow_water_test_factory import TestConservativeGradientJump2D3NElement
from shallow_water_test_factory import TestConservativeFluxCorrected2D3NElement
from shallow_water_test_factory import TestPrimitive2D3NElement
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
from shallow_water_test_factory import TestDamBreakValidation
from shallow_water_test_factory import TestMacDonaldShockValidation
from shallow_water_test_factory import TestSolitaryWaveValidation
from processes_tests.test_line_graph_output_process import TestLineGraphOutputProcess
from processes_tests.test_derivatives_recovery_process import TestDerivativesRecoveryProcess
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
    smallSuite.addTests(TestLoader().loadTestsFromTestCase(TestConservativeResidualViscosity2D3NElement))
    smallSuite.addTests(TestLoader().loadTestsFromTestCase(TestConservativeGradientJump2D3NElement))
    smallSuite.addTests(TestLoader().loadTestsFromTestCase(TestConservativeFluxCorrected2D3NElement))
    smallSuite.addTests(TestLoader().loadTestsFromTestCase(TestPrimitive2D3NElement))
    # smallSuite.addTests(TestLoader().loadTestsFromTestCase(TestBoussinesq2D3NElement)) # FIXME: This test is failing randomly with clang
    smallSuite.addTests(TestLoader().loadTestsFromTestCase(TestSetTopographyProcess))
    smallSuite.addTests(TestLoader().loadTestsFromTestCase(TestVisualizationMeshProcess))
    smallSuite.addTests(TestLoader().loadTestsFromTestCase(TestMacDonaldShockBenchmark))
    smallSuite.addTests(TestLoader().loadTestsFromTestCase(TestMacDonaldTransitionBenchmark))
    smallSuite.addTests(TestLoader().loadTestsFromTestCase(TestDamBreakBenchmark))
    smallSuite.addTests(TestLoader().loadTestsFromTestCase(TestDryDamBreakBenchmark))
    smallSuite.addTests(TestLoader().loadTestsFromTestCase(TestPlanarSurfaceInParabolaBenchmark))
    smallSuite.addTests(TestLoader().loadTestsFromTestCase(TestSolitaryWaveBenchmark))
    smallSuite.addTests(TestLoader().loadTestsFromTestCase(TestLineGraphOutputProcess))
    smallSuite.addTests(TestLoader().loadTestsFromTestCase(TestDerivativesRecoveryProcess))
    smallSuite.addTests(TestLoader().loadTestsFromTestCase(TestWaveGeneratorProcess))

    # Create a test suit with the selected tests plus all small tests
    nightlySuite = suites['nightly']
    nightlySuite.addTests(smallSuite)
    nightlySuite.addTests(TestLoader().loadTestsFromTestCase(TestMeshMovingStrategy))

    # Create a test suit with the validation tests plus all the nightly tests
    validationSuite = suites['validation']
    validationSuite.addTests(nightlySuite)
    validationSuite.addTests(TestLoader().loadTestsFromTestCase(TestDamBreakValidation))
    validationSuite.addTests(TestLoader().loadTestsFromTestCase(TestMacDonaldShockValidation))
    validationSuite.addTests(TestLoader().loadTestsFromTestCase(TestSolitaryWaveValidation))

    # Create a test suit that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(validationSuite)

    return suites

def run():
    KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
    KratosUnittest.runTests(AssembleTestSuites())

if __name__ == '__main__':
    run()
