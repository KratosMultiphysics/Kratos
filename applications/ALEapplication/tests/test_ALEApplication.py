# import Kratos
from KratosMultiphysics import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.ALEApplication import *

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
## SMALL TESTS
from SmallTests import MeshTest2D3NLaplacianTest as TMeshTest2D3NLaplacianTest
from SmallTests import MeshTest2D4NLaplacianTest as TMeshTest2D4NLaplacianTest
from SmallTests import MeshTest3D4NLaplacianTest as TMeshTest3D4NLaplacianTest
from SmallTests import MeshTest2D3NStructuralSimilarityTest as TMeshTest2D3NStructuralSimilarityTest
from SmallTests import MeshTest2D4NStructuralSimilarityTest as TMeshTest2D4NStructuralSimilarityTest
from SmallTests import MeshTest3D4NStructuralSimilarityTest as TMeshTest3D4NStructuralSimilarityTest

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
    # smallSuite.addTest(TMeshTest2D3NLaplacianTest('test_execution'))
    # smallSuite.addTest(TMeshTest2D4NLaplacianTest('test_execution'))
    smallSuite.addTest(TMeshTest2D3NStructuralSimilarityTest('test_execution'))
    smallSuite.addTest(TMeshTest2D4NStructuralSimilarityTest('test_execution'))

    # Create a test suit with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)

    # For very long tests that should not be in nighly and you can use to validate
    validationSuite = suites['validation']
    validationSuite.addTests(smallSuite)
    # smallSuite.addTest(TMeshTest3D4NLaplacianTest('test_execution'))
    smallSuite.addTest(TMeshTest3D4NStructuralSimilarityTest('test_execution'))

    # Create a test suit that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            # TMeshTest2D3NLaplacianTest,
            # TMeshTest2D4NLaplacianTest,
            TMeshTest2D3NStructuralSimilarityTest,
            TMeshTest2D4NStructuralSimilarityTest,
        ])
    )

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssambleTestSuites())
