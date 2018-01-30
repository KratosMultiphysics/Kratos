# import Kratos
import KratosMultiphysics
import KratosMultiphysics.MeshingApplication         as MeshingApplication
try:
  import KratosMultiphysics.ExternalSolversApplication as ExternalSolversApplication
  missing_external_solver_dependencies = False
except ImportError as e:
    missing_external_solver_dependencies = True
try:
  import KratosMultiphysics.FluidDynamicsApplication as FluidDynamicsApplication
  missing_external_fluid_dependencies = False
except ImportError as e:
    missing_external_fluid_dependencies = True
try:
  import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
  missing_external_solid_dependencies = False
except ImportError as e:
    missing_external_solid_dependencies = True

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
## SMALL TESTS 
from test_refine import TestRedistance                      as TTestRedistance
from test_remesh_sphere import TestRemeshMMG                as TTestRemeshMMG
from SmallTests  import TwoDHessianTest                     as TTwoDHessianTest
from SmallTests  import ThreeDHessianTest                   as TThreeDHessianTest
from SmallTests  import TwoDCavityTest                      as TTwoDCavityTest
from SmallTests  import TwoDDynamicBeamTest                 as TTwoDDynamicBeamTest
from SmallTests  import TwoDDynamicBeamLineLoadTest         as TTwoDDynamicBeamLineLoadTest
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
    if( hasattr(MeshingApplication,  "TetrahedraReconnectUtility") ):
        smallSuite.addTest(TTestRedistance('test_refine_all'))
        smallSuite.addTest(TTestRedistance('test_refine_half'))
        smallSuite.addTest(TTestRedistance('test_refine_half_and_improve'))
    else:
        print("TetrahedraReconnectUtility process is not compiled and the corresponding tests will not be executed")
    if( hasattr(MeshingApplication,  "MmgProcess2D") ):
        if (missing_external_fluid_dependencies == False):
            smallSuite.addTest(TTwoDHessianTest('test_execution'))
            smallSuite.addTest(TThreeDHessianTest('test_execution'))
            smallSuite.addTest(TTwoDCavityTest('test_execution'))
            smallSuite.addTest(TTestRemeshMMG('test_remesh_sphere'))
        if (missing_external_solid_dependencies == False):
            smallSuite.addTest(TTwoDDynamicBeamTest('test_execution'))
            smallSuite.addTest(TTwoDDynamicBeamLineLoadTest('test_execution'))
            smallSuite.addTest(TThreeDDynamicBeamTest('test_execution'))
            smallSuite.addTest(TTwoDDynamicPlasticBeamTest('test_execution'))
    else:
        print("MMG process is not compiled and the corresponding tests will not be executed")

    # Create a test suit with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)
    if( hasattr(MeshingApplication,  "MmgProcess2D") ):
        if (missing_external_fluid_dependencies == False):
            nightSuite.addTest(TStanfordBunnyTest('test_execution'))
    else:
        print("MMG process is not compiled and the corresponding tests will not be executed")
    
    # For very long tests that should not be in nighly and you can use to validate 
    validationSuite = suites['validation']
    if( hasattr(MeshingApplication,  "MmgProcess2D") ):
        if (missing_external_fluid_dependencies == False):
            validationSuite.addTest(TTwoDSphereRemeshedChannelTest('test_execution'))
            validationSuite.addTest(TThreeDSphereRemeshedChannelTest('test_execution'))
    else:
        print("MMG process is not compiled and the corresponding tests will not be executed")

    # Create a test suit that contains all the tests:
    allSuite = suites['all']
    if( hasattr(MeshingApplication,  "TetrahedraReconnectUtility") ):
        allSuite.addTests(
            KratosUnittest.TestLoader().loadTestsFromTestCases([
                TTestRedistance
            ])
        )
    else:
        print("TetrahedraReconnectUtility process is not compiled and the corresponding tests will not be executed")
        
    if( hasattr(MeshingApplication,  "MmgProcess2D") ):
        if (missing_external_fluid_dependencies == False):
            allSuite.addTests(
                KratosUnittest.TestLoader().loadTestsFromTestCases([
                    TTwoDHessianTest,
                    TThreeDHessianTest,
                    TTwoDCavityTest,
                    TTestRemeshMMG,
                    #TStanfordBunnyTest,
                    #TTwoDSphereRemeshedChannelTest,
                    #TThreeDSphereRemeshedChannelTest,
                ])
            )
        if (missing_external_solid_dependencies == False):
            allSuite.addTests(
                KratosUnittest.TestLoader().loadTestsFromTestCases([
                    TTwoDDynamicBeamTest,
                    TTwoDDynamicBeamLineLoadTest,
                    TThreeDDynamicBeamTest,
                    ##TTwoDDynamicPlasticBeamTest,
                ])
            )
    else:
        print("MMG process is not compiled and the corresponding tests will not be executed")

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssambleTestSuites())
