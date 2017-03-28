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
from TestExactIntegration import TestLineExactIntegration1 as TTestLineExactIntegration1
from TestExactIntegration import TestLineExactIntegration2 as TTestLineExactIntegration2
from TestExactIntegration import TestTriangleExactIntegration1 as TTestTriangleExactIntegration1
from TestExactIntegration import TestTriangleExactIntegration2 as TTestTriangleExactIntegration2
from TestExactIntegration import TestTriangleExactIntegration3 as TTestTriangleExactIntegration3
from TestExactIntegration import TestQuadrilateralExactIntegration1 as TTestQuadrilateralExactIntegration1
from TestExactIntegration import TestQuadrilateralExactIntegration2 as TTestQuadrilateralExactIntegration2
from SmallTests import SimplePatchTestTwoDMeshTying as TSimplePatchTestTwoDMeshTying
from SmallTests import SimplePatchTestThreeDMeshTying as TSimplePatchTestThreeDMeshTying
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
from SmallTests import ALMHyperSimplePatchTestContact as TALMHyperSimplePatchTestContact
from SmallTests import ALMSimplePatchTestContact as TALMSimplePatchTestContact
from SmallTests import ALMSimplestPatchTestThreeDContact as TALMSimplestPatchTestThreeDContact
from SmallTests import ALMSimplePatchTestThreeDContact as TALMSimplePatchTestThreeDContact

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
    smallSuite.addTest(TTestLineExactIntegration1('test_execution'))
    smallSuite.addTest(TTestLineExactIntegration2('test_execution'))
    smallSuite.addTest(TTestTriangleExactIntegration1('test_execution'))
    smallSuite.addTest(TTestTriangleExactIntegration2('test_execution'))
    smallSuite.addTest(TTestTriangleExactIntegration3('test_execution'))
    smallSuite.addTest(TTestQuadrilateralExactIntegration1('test_execution'))
    smallSuite.addTest(TTestQuadrilateralExactIntegration2('test_execution'))
    smallSuite.addTest(TSimplePatchTestTwoDMeshTying('test_execution'))
    smallSuite.addTest(TSimplePatchTestThreeDMeshTying('test_execution'))
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
    smallSuite.addTest(TALMHyperSimplePatchTestContact('test_execution'))
    smallSuite.addTest(TALMSimplePatchTestContact('test_execution'))
    smallSuite.addTest(TALMSimplestPatchTestThreeDContact('test_execution'))
    smallSuite.addTest(TALMSimplePatchTestThreeDContact('test_execution'))

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
            TTestLineExactIntegration1,
            TTestLineExactIntegration2,
            TTestTriangleExactIntegration1,
            TTestTriangleExactIntegration2,
            TTestTriangleExactIntegration3,
            TTestQuadrilateralExactIntegration1,
            TTestQuadrilateralExactIntegration2,
            TSimplePatchTestTwoDMeshTying,
            #TSimplePatchTestThreeDMeshTying,
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
            TALMHyperSimplePatchTestContact,
            TALMSimplePatchTestContact,
            TALMSimplestPatchTestThreeDContact,
            TALMSimplePatchTestThreeDContact,
            ############# JUST TESTING ###########
            ##THertzCompleteTestContact,
            ##TIroningTestContact,
            ##TIroningDieTestContact,
        ])
    )

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssambleTestSuites())
