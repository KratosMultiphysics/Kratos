# import Kratos
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication
import KratosMultiphysics.RANSModellingApplication

import subprocess

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suites
from evm_k_epsilon_tests import EvmKEpsilonTest
from custom_process_tests import CustomProcessTest

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

    # Create a test suite with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)
    nightSuite.addTest(CustomProcessTest('testApplyFlagProcess'))
    nightSuite.addTest(CustomProcessTest('testScalarCellCenterAveragingProcess'))
    nightSuite.addTest(CustomProcessTest('testVectorCellCenterAveragingProcess'))
    nightSuite.addTest(CustomProcessTest('testVectorAlignProcessTangential'))
    nightSuite.addTest(CustomProcessTest('testVectorAlignProcessNormal'))
    nightSuite.addTest(CustomProcessTest('testWallDistanceCalculationProcess'))
    nightSuite.addTest(CustomProcessTest('testLogarithmicYPlusCalculationProcess'))
    nightSuite.addTest(CustomProcessTest('testNutKEpsilonHighReCalculationProcess'))
    nightSuite.addTest(EvmKEpsilonTest('testBackwardFacingStepKEpsilonTransient'))

    # nightSuite.addTest(EvmKEpsilonTest('testCylinderTransient')) #TODO: Has a convergence problem, therefore gives a race condition

    # For very long tests that should not be in nighly and you can use to validate
    validationSuite = suites['validation']
    validationSuite.addTest(EvmKEpsilonTest('testChannelTransient'))

    # Create a test suite that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(nightSuite)

    # allSuite.addTests(validationSuite)

    return suites

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning cpp unit tests ...")
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished running cpp unit tests!")

    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning mpi python tests ...")
    try:
        import KratosMultiphysics.mpi as KratosMPI
        import KratosMultiphysics.MetisApplication as MetisApplication
        import KratosMultiphysics.TrilinosApplication as TrilinosApplication
        p = subprocess.Popen(["mpiexec", "-np", "2", "python3", "test_FluidDynamicsApplication_mpi.py"], stdout=subprocess.PIPE)
        p.wait()
        KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished mpi python tests!")
    except ImportError:
        KratosMultiphysics.Logger.PrintInfo("Unittests", "mpi is not available!")

    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning python tests ...")
    KratosUnittest.runTests(AssembleTestSuites())
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished python tests!")