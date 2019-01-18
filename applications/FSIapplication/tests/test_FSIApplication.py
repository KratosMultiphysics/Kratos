# import Kratos
from KratosMultiphysics import *
from KratosMultiphysics.FSIApplication import *

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
## SMALL TESTS
from convergence_accelerator_test import ConvergenceAcceleratorTest
from convergence_accelerator_spring_test import ConvergenceAcceleratorSpringTest
from FSI_problem_emulator_test import FSIProblemEmulatorTest
from non_conformant_one_side_map_test import NonConformantOneSideMapTest

## NIGTHLY TESTS

## VALIDATION TESTS
from mok_benchmark_test import MokBenchmarkTest

def AssembleTestSuites():
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
    smallSuite.addTest(ConvergenceAcceleratorTest('test_aitken_accelerator'))
    smallSuite.addTest(ConvergenceAcceleratorTest('test_mvqn_accelerator'))
    smallSuite.addTest(ConvergenceAcceleratorTest('test_mvqn_recusive_accelerator'))
    smallSuite.addTest(ConvergenceAcceleratorTest('test_accelerator_with_jacobian'))
    smallSuite.addTest(FSIProblemEmulatorTest('testFSIProblemEmulatorWithAitken'))
    smallSuite.addTest(FSIProblemEmulatorTest('testFSIProblemEmulatorWithMVQN'))
    smallSuite.addTest(FSIProblemEmulatorTest('testFSIProblemEmulatorWithMVQNRecursive'))
    smallSuite.addTest(ConvergenceAcceleratorSpringTest('test_aitken_accelerator_constant_forces'))
    smallSuite.addTest(ConvergenceAcceleratorSpringTest('test_aitken_accelerator_variable_stiffness'))
    smallSuite.addTest(ConvergenceAcceleratorSpringTest('test_mvqn_recursive_accelerator_constant_forces'))
    smallSuite.addTest(ConvergenceAcceleratorSpringTest('test_mvqn_recursive_accelerator_variable_stiffness'))
    smallSuite.addTest(NonConformantOneSideMapTest('test2D_1'))
    smallSuite.addTest(NonConformantOneSideMapTest('test2D_2'))
    smallSuite.addTest(NonConformantOneSideMapTest('test3D_1'))
    smallSuite.addTest(NonConformantOneSideMapTest('test3D_two_faces'))

    # Create a test suit with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)

    # For very long tests that should not be in nighly and you can use to validate
    validationSuite = suites['validation']
    validationSuite.addTest(MokBenchmarkTest('testMokBenchmark'))

    # Create a test suit that contains all the tests
    allSuite = suites['all']
    allSuite.addTests(nightSuite)

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
