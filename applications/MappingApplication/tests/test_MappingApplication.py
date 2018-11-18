# import Kratos
import KratosMultiphysics
import KratosMultiphysics.MappingApplication

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

import run_cpp_unit_tests

# Import the tests or test_classes to create the suits
from SmallTests import NearestNeighborTest_1 as TNearestNeighborTest_1
from SmallTests import NearestElementTest2D_1 as TNearestElementTest2D_1
from SmallTests import MapperTests as TMapperTests
from test_mapper_tests import MapperTests
from test_mapper_flags_tests import MapperFlagsTests

KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
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
    smallSuite = suites['small']
    smallSuite.addTest(TNearestNeighborTest_1('test_execution'))
    smallSuite.addTest(TNearestElementTest2D_1('test_execution'))
    smallSuite.addTest(TMapperTests('test_execution'))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([MapperTests]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([MapperFlagsTests]))


    # Create a test suit with the selected tests
    # nightSuite will contain the following tests:
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)

    # Create a test suit that contains all the tests from every testCase
    # in the list:
    allSuite = suites['all']
    allSuite.addTests(nightSuite) # already contains the smallSuite

    return suites

if __name__ == '__main__':
    print("Running cpp unit tests for MappingApplication...")
    run_cpp_unit_tests.run()
    print("Finished running cpp unit tests!")
    print("Running python tests for MappingApplication...")
    KratosUnittest.runTests(AssembleTestSuites())
    print("Finished python tests!")