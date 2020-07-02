import subprocess
import os.path

# import Kratos
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication
import KratosMultiphysics.RANSApplication

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

import KratosMultiphysics.kratos_utilities as kratos_utilities

# Import the tests o test_classes to create the suites
from evm_k_epsilon_tests import EvmKEpsilonTest
from custom_process_tests import CustomProcessTest
from adjoint_k_epsilon_sensitivity_2d import AdjointKEpsilonSensitivity2D
import run_cpp_unit_tests


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

    # Create a test suite with the selected tests (Small tests):
    smallSuite = suites['small']
    smallSuite.addTest(EvmKEpsilonTest('testOneElementKEpsilonSteady'))
    smallSuite.addTest(AdjointKEpsilonSensitivity2D('testOneElementSteady'))

    # Create a test suite with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)

    nightSuite.addTest(CustomProcessTest('testCheckScalarBoundsProcess'))
    nightSuite.addTest(CustomProcessTest('testCheckVectorBoundsProcess'))
    nightSuite.addTest(CustomProcessTest('testClipScalarVariableProcess'))
    nightSuite.addTest(CustomProcessTest('testApplyFlagProcess'))
    nightSuite.addTest(CustomProcessTest('testWallDistanceCalculationProcess'))
    nightSuite.addTest(CustomProcessTest('testLogarithmicYPlusCalculationProcess'))
    nightSuite.addTest(CustomProcessTest('testLogarithmicYPlusVelocitySensitivitiesProcessFlow'))
    nightSuite.addTest(CustomProcessTest('testNutKEpsilonHighReCalculationProcess'))
    nightSuite.addTest(CustomProcessTest('testNutKEpsilonHighReSensitivitiesProcess'))
    nightSuite.addTest(EvmKEpsilonTest('testBackwardFacingStepKEpsilonTransient'))
    nightSuite.addTest(EvmKEpsilonTest('testChannelFlowKEpsilonSteady'))
    nightSuite.addTest(EvmKEpsilonTest('testChannelFlowKEpsilonSteadyPeriodic'))

    # For very long tests that should not be in nighly and you can use to validate
    # validationSuite = suites['validation']

    # Create a test suite that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(nightSuite)

    # allSuite.addTests(validationSuite)

    return suites


if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(
        KratosMultiphysics.Logger.Severity.WARNING)
    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning cpp unit tests ...")
    run_cpp_unit_tests.run()
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished running cpp unit tests!")

    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning python tests ...")
    KratosUnittest.runTests(AssembleTestSuites())
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished python tests!")
