# import Kratos
from KratosMultiphysics import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

## SMALL TESTS
#from SmallTests import 

## NIGHTLY TESTS

## VALIDATION TESTS

def AssambleTestSuites():
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
    #smallSuite.addTest()

    # Create a test suit with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)

    # For very long tests that should not be in nightly and you can use to validate
    validationSuite = suites['validation']
    validationSuite.addTests(nightSuite)

    # Create a test suit that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(nightSuite)

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssambleTestSuites())
