# import Kratos
from KratosMultiphysics import *
from KratosMultiphysics.ExternalSolversApplication import *
import run_cpp_unit_tests

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

## SMALL TESTS
from SmallTests import Pfem2PrimitiveVariables as TPfem2PrimitiveVariables

## NIGHTLY TESTS
from NightlyTests import Pfem2ConservedVariables as TPfem2ConservedVariables
from NightlyTests import EulerianPrimitiveVariables as TEulerianPrimitiveVariables
from NightlyTests import EulerianConservedVariables as TEulerianConservedVariables

## VALIDATION TESTS

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
    smallSuite.addTest(TPfem2PrimitiveVariables('test_execution'))

    # Create a test suit with the selected tests plus all small tests
    nightlySuite = suites['nightly']
    nightlySuite.addTests(smallSuite)
    nightlySuite.addTest(TPfem2ConservedVariables('test_execution'))
    # nightlySuite.addTest(TEulerianPrimitiveVariables('test_execution'))
    # nightlySuite.addTest(TEulerianConservedVariables('test_execution'))

    # Create a test suit that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(nightlySuite)

    return suites

if __name__ == '__main__':
    run_cpp_unit_tests.run()
    KratosUnittest.runTests(AssembleTestSuites())
