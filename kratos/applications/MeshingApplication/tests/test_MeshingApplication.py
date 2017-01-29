# import Kratos
import KratosMultiphysics
import KratosMultiphysics.ExternalSolversApplication
import KratosMultiphysics.FluidDynamicsApplication
import KratosMultiphysics.MeshingApplication

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
## SMALL TESTS # TODO: Add test_refine.py to the list of tests
from SmallTests import CoarseSphereTest as TCoarseSphereTest

## SMALL TESTS
from NightlyTests import StanfordBunnyTest as TStanfordBunnyTest

## VALIDATION TESTS 

def AssambleTestSuites():
    ''' Populates the test suites to run.

    Populates the test suites to run. At least, it should pupulate the suites:
    "small", "nighlty" and "all"

    Return
    ------

    suites: A dictionary of suites
        The set of suites with its test_cases added.
    '''
    suites = KratosUnittest.KratosSuites

    # Create a test suit with the selected tests (Small tests):
    smallSuite = suites['small']
    smallSuite.addTest(TCoarseSphereTest('test_execution'))

    # Create a test suit with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)
    nightSuite.addTest(TStanfordBunnyTest('test_execution'))
    
    # For very long tests that should not be in nighly and you can use to validate 
    validationSuite = suites['validation']

    # Create a test suit that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
        ])
    )

    if( hasattr(KratosMultiphysics.MeshingApplication,  "MmgUtility2D") ):
        allSuite.addTests(
            KratosUnittest.TestLoader().loadTestsFromTestCases([
                TCoarseSphereTest,
                TStanfordBunnyTest,
            ])
        )
    else:
        print("MMG utility is not compiled and the corresponding tests will not be executed")

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssambleTestSuites())
