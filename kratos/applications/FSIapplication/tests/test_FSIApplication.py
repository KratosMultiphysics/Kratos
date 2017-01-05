# import Kratos
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.FSIApplication import *

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
## SMALL TESTS
from SmallTests import FSIProblemEmulatorTest as TFSIProblemEmulatorTest
from SmallTests import NonConformantOneSideMap2D_test1 as TNonConformantOneSideMap2D_test1
from SmallTests import NonConformantOneSideMap2D_test2 as TNonConformantOneSideMap2D_test2
from KratosExecuteConvergenceAcceleratorTest import KratosExecuteConvergenceAcceleratorTest as TConvergenceAcceleratorTest
#~ from SmallTests import NonConformantOneSideMap3D_test2 as TNonConformantOneSideMap3D_test2

## NIGTHLY TESTS
#~ from NightlyTests import MokBenchmarkTest as TMokBenchmarkTest

## VALIDATION TESTS
from ValidationTests import MokBenchmarkTest as TMokBenchmark

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
    smallSuite.addTest(TNonConformantOneSideMap2D_test1('test_execution'))
    smallSuite.addTest(TNonConformantOneSideMap2D_test2('test_execution'))
    smallSuite.addTest(TConvergenceAcceleratorTest('test_aitken_accelerator'))
    smallSuite.addTest(TConvergenceAcceleratorTest('test_mvqn_accelerator'))
    smallSuite.addTest(TConvergenceAcceleratorTest('test_mvqn_recusive_accelerator'))
    smallSuite.addTest(TConvergenceAcceleratorTest('test_accelerator_with_jacobian'))
    smallSuite.addTest(TFSIProblemEmulatorTest('test_execution'))
    #~ smallSuite.addTest(TNonConformantOneSideMap3D_test2('test_execution'))

    # Create a test suit with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)
    #~ nightSuite.addTest(TMokBenchmarkTest('test_execution'))

    # For very long tests that should not be in nighly and you can use to validate
    validationSuite = suites['validation']
    validationSuite.addTest(TMokBenchmark('test_execution'))

    # Create a test suit that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            TNonConformantOneSideMap2D_test1,
            TNonConformantOneSideMap2D_test2,
            TConvergenceAcceleratorTest,
            TFSIProblemEmulatorTest,
            TMokBenchmark
            #####TTurekBenchmarkTest
        ])
    )

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssambleTestSuites())
