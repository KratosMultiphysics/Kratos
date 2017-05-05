# import Kratos
import KratosMultiphysics
import KratosMultiphysics.ExternalSolversApplication as ExternalSolversApplication
import KratosMultiphysics.FluidDynamicsApplication   as FluidDynamicsApplication
import KratosMultiphysics.MeshingApplication         as MeshingApplication

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
## SMALL TESTS 
from test_refine import TestRedistance                      as TTestRedistance
from SmallTests  import TwoDHessianTest                     as TTwoDHessianTest
from SmallTests  import ThreeDHessianTest                   as TThreeDHessianTest
from SmallTests  import TwoDCavityTest                      as TTwoDCavityTest
from SmallTests  import CoarseSphereTest                    as TCoarseSphereTest
from SmallTests  import TwoDDynamicBeamTest                 as TTwoDDynamicBeamTest
from SmallTests  import ThreeDDynamicBeamTest               as TThreeDDynamicBeamTest
from SmallTests  import TwoDDynamicPlasticBeamTest          as TTwoDDynamicPlasticBeamTest

## NIGHTLY TESTS
from NightlyTests import StanfordBunnyTest                  as TStanfordBunnyTest

## VALIDATION TESTS 
from ValidationTests import TwoDSphereRemeshedChannelTest   as TTwoDSphereRemeshedChannelTest
from ValidationTests import ThreeDSphereRemeshedChannelTest as TThreeDSphereRemeshedChannelTest

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
    smallSuite.addTest(TTestRedistance('test_refine_all'))
    smallSuite.addTest(TTestRedistance('test_refine_half'))
    smallSuite.addTest(TTestRedistance('test_refine_half_and_improve'))
    if( hasattr(MeshingApplication,  "MmgUtility2D") ):
        smallSuite.addTest(TTwoDHessianTest('test_execution'))
        smallSuite.addTest(TThreeDHessianTest('test_execution'))
        smallSuite.addTest(TTwoDCavityTest('test_execution'))
        smallSuite.addTest(TCoarseSphereTest('test_execution'))
        smallSuite.addTest(TTwoDDynamicBeamTest('test_execution'))
        smallSuite.addTest(TThreeDDynamicBeamTest('test_execution'))
        smallSuite.addTest(TTwoDDynamicPlasticBeamTest('test_execution'))
    else:
        print("MMG utility is not compiled and the corresponding tests will not be executed")

    # Create a test suit with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)
    if( hasattr(MeshingApplication,  "MmgUtility2D") ):
        nightSuite.addTest(TStanfordBunnyTest('test_execution'))
    else:
        print("MMG utility is not compiled and the corresponding tests will not be executed")
    
    # For very long tests that should not be in nighly and you can use to validate 
    validationSuite = suites['validation']
    if( hasattr(MeshingApplication,  "MmgUtility2D") ):
        validationSuite.addTest(TTwoDSphereRemeshedChannelTest('test_execution'))
        validationSuite.addTest(TThreeDSphereRemeshedChannelTest('test_execution'))
    else:
        print("MMG utility is not compiled and the corresponding tests will not be executed")

    # Create a test suit that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            TTestRedistance
        ])
    )

    if( hasattr(MeshingApplication,  "MmgUtility2D") ):
        allSuite.addTests(
            KratosUnittest.TestLoader().loadTestsFromTestCases([
                TTwoDHessianTest,
                TThreeDHessianTest,
                TTwoDCavityTest,
                TCoarseSphereTest,
                TTwoDDynamicBeamTest,
                #TThreeDDynamicBeamTest,
                #TTwoDDynamicPlasticBeamTest,
                #TStanfordBunnyTest,
                #TTwoDSphereRemeshedChannelTest,
                #TThreeDSphereRemeshedChannelTest,
            ])
        )
    else:
        print("MMG utility is not compiled and the corresponding tests will not be executed")

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssambleTestSuites())
