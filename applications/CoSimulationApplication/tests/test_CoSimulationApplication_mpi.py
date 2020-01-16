# Importing the Kratos Library
import KratosMultiphysics as KM

if not KM.IsDistributedRun():
    raise Exception("This test script can only be executed in MPI!")

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import tests
from test_convergence_criteria import TestConvergenceCriteriaWrapper
from test_convergence_accelerators import TestConvergenceAcceleratorWrapper


def AssembleTestSuites():
    ''' Populates the test suites to run.

    Populates the test suites to run. At least, it should populate the suites:
    "mpi_small", "mpi_nighlty" and "mpi_all"

    Return
    ------

    suites: A dictionary of suites
        The set of suites with its test_cases added.
    '''
    suites = KratosUnittest.KratosSuites
    ################################################################################
    smallSuite = suites['mpi_small'] # These tests are executed by the continuous integration tool
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestConvergenceCriteriaWrapper]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestConvergenceAcceleratorWrapper]))

    ################################################################################
    nightSuite = suites['mpi_nightly'] # These tests are executed in the nightly build
    nightSuite.addTests(smallSuite)

    ################################################################################
    # For very long tests that should not be in nighly and you can use to validate
    validationSuite = suites['mpi_validation']

    ################################################################################
    # Create a test suit that contains all the tests:
    allSuite = suites['mpi_all']
    allSuite.addTests(nightSuite) # already contains the smallSuite
    allSuite.addTests(validationSuite)

    return suites


if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
