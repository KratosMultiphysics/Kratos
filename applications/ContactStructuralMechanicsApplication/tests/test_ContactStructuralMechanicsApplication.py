# import Kratos
from KratosMultiphysics import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *
from KratosMultiphysics.ContactStructuralMechanicsApplication import *

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
## SMALL TESTS
#from SmallTests import BasicCATest as TBasicCATest
#from SmallTests import SolidCATest as TSolidCATest
from SmallTests import SimplePatchTestContact as TSimplePatchTestContact
from SmallTests import SimpleSlopePatchTestContact as TSimpleSlopePatchTestContact
from SmallTests import SimplePatchNotMatchingATestContact as TSimplePatchNotMatchingATestContact
from SmallTests import SimplePatchNotMatchingBTestContact as TSimplePatchNotMatchingBTestContact
from SmallTests import TaylorPatchTestContact as TTaylorPatchTestContact
from SmallTests import TaylorPatchDynamicTestContact as TTaylorPatchDynamicTestContact
from SmallTests import HertzSimpleTestContact as THertzSimpleTestContact
from SmallTests import HertzSimpleSphereTestContact as THertzSimpleSphereTestContact
from SmallTests import HertzSphereTestContact as THertzSphereTestContact
#from SmallTests import HertzCompleteTestContact as THertzCompleteTestContact
from SmallTests import ThreeDSimplestPatchMatchingTestContact as TThreeDSimplestPatchMatchingTestContact
from SmallTests import ThreeDSimplestTrianglePatchMatchingTestContact as TThreeDSimplestTrianglePatchMatchingTestContact
from SmallTests import ThreeDPatchMatchingTestContact as TThreeDPatchMatchingTestContact
from SmallTests import ThreeDPatchNotMatchingTestContact as TThreeDPatchNonMatchingTestContact

## NIGTHLY TESTS
from NightlyTests import IroningTestContact as TIroningTestContact
from NightlyTests import IroningDieTestContact as TIroningDieTestContact

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
    #smallSuite.addTest(TBasicCATest('test_execution'))
    #smallSuite.addTest(TSolidCATest('test_execution'))
    smallSuite.addTest(TSimplePatchTestContact('test_execution'))
    smallSuite.addTest(TSimpleSlopePatchTestContact('test_execution'))
    smallSuite.addTest(TSimplePatchNotMatchingATestContact('test_execution'))
    smallSuite.addTest(TSimplePatchNotMatchingBTestContact('test_execution'))
    smallSuite.addTest(TTaylorPatchTestContact('test_execution'))
    smallSuite.addTest(TTaylorPatchDynamicTestContact('test_execution'))
    smallSuite.addTest(THertzSimpleSphereTestContact('test_execution'))
    smallSuite.addTest(THertzSphereTestContact('test_execution'))
    smallSuite.addTest(THertzSimpleTestContact('test_execution'))
    #smallSuite.addTest(THertzCompleteTestContact('test_execution'))
    smallSuite.addTest(TThreeDSimplestPatchMatchingTestContact('test_execution'))
    smallSuite.addTest(TThreeDSimplestTrianglePatchMatchingTestContact('test_execution'))
    smallSuite.addTest(TThreeDPatchMatchingTestContact('test_execution'))
    smallSuite.addTest(TThreeDPatchNonMatchingTestContact('test_execution'))

    # Create a test suit with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)
    #nightSuite.addTest(TIroningTestContact('test_execution'))
    #nightSuite.addTest(TIroningDieTestContact('test_execution'))
    
    # For very long tests that should not be in nighly and you can use to validate 
    validationSuite = suites['validation']

    # Create a test suit that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            #TBasicCATest,
            #TSolidCATest,
            TSimplePatchTestContact,
            TSimpleSlopePatchTestContact,
            TSimplePatchNotMatchingATestContact,
            TSimplePatchNotMatchingBTestContact,
            TTaylorPatchTestContact,
            TTaylorPatchDynamicTestContact,
            THertzSimpleTestContact,
            THertzSimpleSphereTestContact,
            THertzSphereTestContact,
            TThreeDSimplestPatchMatchingTestContact,
            TThreeDSimplestTrianglePatchMatchingTestContact,
            TThreeDPatchMatchingTestContact,
            TThreeDPatchNonMatchingTestContact,
            ############# JUST TESTING ###########
            ##THertzCompleteTestContact,
            ##TIroningTestContact,
            ##TIroningDieTestContact,
        ])
    )

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssambleTestSuites())
