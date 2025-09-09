# import Kratos
import KratosMultiphysics
import KratosMultiphysics.MetisApplication

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
from test_quad_partition import TestQuadPartition

def AssembleTestSuites():
    ''' Populates the test suites to run.

    Populates the test suites to run. At least, it should pupulate the suites:
    "small", "nightly" and "all"

    Return
    ------

    suites: A dictionary of suites
        The set of suites with its test_cases added.
    '''

    suites = KratosUnittest.KratosSuites

    # Create a test suit with the selected tests (Small tests):
    # smallSuite will contain the following tests:
    # - testSmallExample
    smallMPISuite = suites['small']
    smallMPISuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestQuadPartition]))

    # Create a test suit with the selected tests
    # nightSuite will contain the following tests:
    # - testSmallExample
    # - testNightlyFirstExample
    # - testNightlySecondExample
    nightMPISuite = suites['nightly']
    nightMPISuite.addTests(smallMPISuite)

    # Create a test suit that contains all the tests from every testCase
    # in the list:
    allMPISuite = suites['all']
    allMPISuite.addTests(nightMPISuite)

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
