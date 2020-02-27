# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
from cable_net_test_factory import TestCableNetCoSimulationCases
from cable_net_test_factory import TestCableNetFEMCases
from test_empirical_spring  import EmpiricalSpringTests
from test_edge_cable_process  import EdgeCableProcessTests

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
    ################################################################################
    smallSuite = suites['small'] # These tests are executed by the continuous integration tool
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCableNetFEMCases]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([EmpiricalSpringTests]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([EdgeCableProcessTests]))


    ################################################################################
    nightSuite = suites['nightly'] # These tests are executed in the nightly build
    nightSuite.addTests(smallSuite)


    ################################################################################
    # For very long tests that should not be in nighly and you can use to validate
    validationSuite = suites['validation']
    validationSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCableNetCoSimulationCases]))


    ################################################################################
    # Create a test suit that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(nightSuite) # already contains the smallSuite
    allSuite.addTests(validationSuite)

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
