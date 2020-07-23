# import Kratos
import KratosMultiphysics
import KratosMultiphysics.ExaquteSandboxApplication

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
from generalTests import KratosExaquteSandboxGeneralTests
from test_divergencefree_refinement_process import TimeAveragingProcessTests

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
    # smallSuite will contain the following tests:
    # - testSmallExample
    smallSuite = suites['small']
    smallSuite.addTest(KratosExaquteSandboxGeneralTests('testSmallExample'))
    smallSuite.addTest(TimeAveragingProcessTests('testDivergenceFreeRefinementProcess'))

    # Create a test suit with the selected tests
    # nightSuite will contain the following tests:
    # - testSmallExample
    # - testNightlyFirstExample
    # - testNightlySecondExample
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)

    # Create a test suit that contains all the tests from every testCase
    # in the list:
    allSuite = suites['all']
    allSuite.addTests(smallSuite)
    allSuite.addTests(nightSuite)

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
