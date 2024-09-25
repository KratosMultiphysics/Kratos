# import Kratos
from KratosMultiphysics import *
from KratosMultiphysics.FSIApplication import *

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
## SMALL TESTS
from convergence_accelerator_test import ConvergenceAcceleratorTest
from convergence_accelerator_spring_test import ConvergenceAcceleratorSpringTest
from embedded_fsi_test import EmbeddedFsiTest
from fsi_coupling_interface_test import FSICouplingInterfaceTest
from FSI_problem_emulator_test import FSIProblemEmulatorTest

## NIGTHLY TESTS

## VALIDATION TESTS
from mok_benchmark_test import MokBenchmarkTest

def AssembleTestSuites():
    ''' Populates the test suites to run.

    Populates the test suites to run. At least, it should populate the suites:
    "small", "nighlty" and "all"

    Return
    ------

    suites: A dictionary of suites
        The set of suites with its test_cases added.
    '''
    suites = KratosUnittest.KratosSuites

    # Create a test suit with the selected tests (Small tests):
    smallSuite = suites['small']
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([ConvergenceAcceleratorTest]))
    # smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([EmbeddedFsiTest])) #TODO: To be activated once #12681 is merged
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([FSICouplingInterfaceTest]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([FSIProblemEmulatorTest]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([ConvergenceAcceleratorSpringTest]))

    # Create a test suit with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)

    # For very long tests that should not be in nighly and you can use to validate
    validationSuite = suites['validation']
    validationSuite.addTest(MokBenchmarkTest('testMokBenchmark'))

    # Create a test suit that contains all the tests
    allSuite = suites['all']
    allSuite.addTests(nightSuite)

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
