# import Kratos
from KratosMultiphysics import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
## SMALL TESTS
from SmallTests import EmbeddedArtificialCompressibilityTest as TEmbeddedArtificialCompressibilityTest
from SmallTests import EmbeddedCouetteTest as TEmbeddedCouetteTest
from SmallTests import EmbeddedCouetteImposedTest as TEmbeddedCouetteImposedTest
from SmallTests import EmbeddedReservoirTest as TEmbeddedReservoirTest
from SmallTests import ManufacturedSolutionTest as TManufacturedSolutionTest
from SmallTests import NavierStokesWallConditionTest as TNavierStokesWallConditionTest

## NIGTHLY TESTS
#~ from NightlyTests import MokBenchmarkTest as TMokBenchmarkTest

## VALIDATION TESTS
#~ from ValidationTests import TurekBenchmarkTest as TTurekBenchmarkTest

def AssambleTestSuites():
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
    smallSuite.addTest(TEmbeddedArtificialCompressibilityTest('test_execution'))
    smallSuite.addTest(TEmbeddedCouetteTest('test_execution'))
    smallSuite.addTest(TEmbeddedCouetteImposedTest('test_execution'))
    smallSuite.addTest(TEmbeddedReservoirTest('test_execution'))
    smallSuite.addTest(TManufacturedSolutionTest('test_execution'))
    smallSuite.addTest(TNavierStokesWallConditionTest('test_execution'))

    # Create a test suit with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)

    # For very long tests that should not be in nighly and you can use to validate
    validationSuite = suites['validation']
    validationSuite.addTests(smallSuite)

    # Create a test suit that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            TEmbeddedArtificialCompressibilityTest,
            TEmbeddedCouetteTest,
            TEmbeddedCouetteImposedTest,
            TEmbeddedReservoirTest,
            TManufacturedSolutionTest,
            TNavierStokesWallConditionTest
        ])
    )

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssambleTestSuites())
