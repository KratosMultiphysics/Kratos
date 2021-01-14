# Import Kratos
import KratosMultiphysics
import KratosMultiphysics.MultilevelMonteCarloApplication
try:
    import KratosMultiphysics.mpi as KratosMPI
except ImportError:
    raise Exception("KratosMPI could not be imported!")

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
    smallMPISuite = suites['mpi_small']

    # Create a test suit with the selected tests plus all small  mpi tests
    nightlyMPISuite = suites['mpi_nightly']
    nightlyMPISuite.addTests(smallMPISuite)
    nightlyMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestXMCAlgorithmMPI]))

    # Create a test suit that contains all the tests
    allMPISuite = suites['mpi_all']
    allMPISuite.addTests(nightlyMPISuite) # already contains the smallMPISuite

    # can be removed after the cmd-line of the testing accepts "--using-mpi"
    allSuite = suites['all']
    allSuite.addTests(allMPISuite)

    return suites

if __name__ == '__main__':
    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning python tests ...")
    KratosUnittest.runTests(AssembleTestSuites())
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished python tests!")
