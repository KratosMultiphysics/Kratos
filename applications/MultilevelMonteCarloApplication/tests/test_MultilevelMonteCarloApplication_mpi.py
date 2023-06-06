# Import Kratos
import KratosMultiphysics
import KratosMultiphysics.MultilevelMonteCarloApplication
if not KratosMultiphysics.IsDistributedRun():
    raise Exception("This test script can only be executed in MPI!")

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import test classes to create the suits
from test_xmcAlgorithm_mpi import TestXMCAlgorithmMPI

def AssembleTestSuites():
    """
    Populates the test suites to run.

    Populates the test suites to run. At least, it should populate the suites:
    "small", "nighlty" and "all".

    Return
    ------

    suites: A dictionary of suites
            The set of suites with its test_cases added.
    """

    suites = KratosUnittest.KratosSuites

    # Create a test suit with the selected tests
    smallSuite = suites['mpi_small']

    # Create a test suit with the selected tests plus all small mpi tests
    nightlySuite = suites['mpi_nightly']
    nightlySuite.addTests(smallSuite)
    nightlySuite.addTest(TestXMCAlgorithmMPI('mpi_test_mlmc_Kratos_ParMmg'))

    # For very long tests that should not be in nightly and you can use to validate
    validationSuite = suites['mpi_validation']
    validationSuite.addTest(TestXMCAlgorithmMPI('mpi_test_mc_Kratos'))
    validationSuite.addTest(TestXMCAlgorithmMPI('mpi_test_mlmc_Kratos'))

    # Create a test suit that contains all the tests
    allSuite = suites['mpi_all']
    allSuite.addTests(nightlySuite) # already contains the smallSuite
    allSuite.addTests(validationSuite)

    return suites

if __name__ == '__main__':
    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning python tests ...")
    KratosUnittest.runTests(AssembleTestSuites())
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished python tests!")
