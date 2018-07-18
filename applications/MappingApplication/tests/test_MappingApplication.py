# import Kratos
from KratosMultiphysics import *

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
from generalTests import KratosMappingGeneralTests

from SmallTests import NearestNeighborTest_1 as TNearestNeighborTest_1
from SmallTests import NearestElementTest2D_1 as TNearestElementTest2D_1
from SmallTests import MapperTests as TMapperTests


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
    smallSuite.addTest(TNearestNeighborTest_1('test_execution'))
    smallSuite.addTest(TNearestElementTest2D_1('test_execution'))
    smallSuite.addTest(TMapperTests('test_execution'))

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
    allSuite.addTests(nightSuite) # already contains the smallSuite

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
