# import Kratos
import KratosMultiphysics
import KratosMultiphysics.MeshingApplication         as MeshingApplication
import run_cpp_unit_tests
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
from SmallTests  import TwoDDynamicBeamTest                 as TTwoDDynamicBeamTest
from SmallTests  import TwoDDynamicBeamLineLoadTest         as TTwoDDynamicBeamLineLoadTest
from SmallTests  import ThreeDDynamicBeamTest               as TThreeDDynamicBeamTest

## NIGHTLY TESTS

## VALIDATION TESTS

def AssembleTestSuites():
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
        smallSuite.addTest(TTestRemeshMMG('test_remesh_sphere'))
        if (missing_external_solid_dependencies is False):
            smallSuite.addTest(TTwoDDynamicBeamTest('test_execution'))
            smallSuite.addTest(TTwoDDynamicBeamLineLoadTest('test_execution'))
            smallSuite.addTest(TThreeDDynamicBeamTest('test_execution'))
    else:
        print("MMG process is not compiled and the corresponding tests will not be executed")

    # Create a test suit with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)
    #if( hasattr(MeshingApplication,  "MmgProcess2D") ):
        #if (missing_external_fluid_dependencies is False):
            #nightSuite.addTest()
    #else:
        #print("MMG process is not compiled and the corresponding tests will not be executed")

    # For very long tests that should not be in nighly and you can use to validate
    validationSuite = suites['validation']
    #if( hasattr(MeshingApplication,  "MmgProcess2D") ):
        #if (missing_external_fluid_dependencies is False):
            #validationSuite.addTest()
    #else:
        #print("MMG process is not compiled and the corresponding tests will not be executed")

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
        allSuite.addTests(
            KratosUnittest.TestLoader().loadTestsFromTestCases([
                TTestRemeshMMG,
            ])
        )
        if (missing_external_solid_dependencies is False):
            allSuite.addTests(
                KratosUnittest.TestLoader().loadTestsFromTestCases([
                    TTwoDDynamicBeamTest,
                    TTwoDDynamicBeamLineLoadTest,
                    TThreeDDynamicBeamTest,
                ])
            )
    else:
        print("MMG process is not compiled and the corresponding tests will not be executed")

    return suites

if __name__ == '__main__':
    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning cpp unit tests ...")
    run_cpp_unit_tests.run()
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished running cpp unit tests!")

    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning python tests ...")
    KratosUnittest.runTests(AssembleTestSuites())
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished python tests!")
