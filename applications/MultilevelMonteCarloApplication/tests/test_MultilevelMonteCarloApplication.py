# import Kratos
import KratosMultiphysics
import KratosMultiphysics.MultilevelMonteCarloApplication as KratosMLMC
import KratosMultiphysics.MeshingApplication as MeshingApplication

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
from test_multilevel_montecarlo import KratosMultilevelMonteCarloGeneralTests


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
    # smallSuite will contain the following tests:
    # -
    smallSuite = suites['small']

    # Create a test suit with the selected tests
    # nightSuite will contain the following tests:
    # - testMonteCarloAnalysis
    # - testMultilevelMonteCarloAnalysis
    nightSuite = suites['nightly']
    if(hasattr(MeshingApplication,"MmgProcess2D")):
        nightSuite.addTest(KratosMultilevelMonteCarloGeneralTests('testMonteCarloAnalysis'))
        nightSuite.addTest(KratosMultilevelMonteCarloGeneralTests('testMultilevelMonteCarloAnalysis'))
    else:
        print("MMG process is not compiled and the corresponding tests will not be executed")

    # Create a test suit that contains all the tests from every testCase
    # in the list:
    allSuite = suites['all']
    allSuite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            KratosMultilevelMonteCarloGeneralTests
        ])
    )

    return suites

if __name__ == '__main__':
    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning python tests ...")
    KratosUnittest.runTests(AssembleTestSuites())
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished python tests!")
