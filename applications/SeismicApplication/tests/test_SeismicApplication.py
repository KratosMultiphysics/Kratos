# import Kratos
import KratosMultiphysics
import KratosMultiphysics.SeismicApplication

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
from generalTests import KratosSeismicGeneralTests
from seismic_test_factory import FiberBeamElementTest


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
    # These tests are executed by the continuous integration tool, so they have to be very fast!
    # Execution time << 1 sec on a regular PC !!!
    # If the tests in the smallSuite take too long then merging to master will not be possible!
    smallSuite = suites['small'] # These tests are executed by the continuous integration tool
    nightSuite = suites['nightly'] # These tests are executed in the nightly build

    # smallSuite.addTest(KratosSeismicGeneralTests('testSmallExample'))
    smallSuite.addTest(FiberBeamElementTest('test_execution'))

    nightSuite.addTests(smallSuite)
    # Create a test suit that contains all the tests from every testCase
    # in the list:
    allSuite = suites['all']
    allSuite.addTests(nightSuite) # Already contains the smallSuite
    # allSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([KratosSeismicGeneralTests]))

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
