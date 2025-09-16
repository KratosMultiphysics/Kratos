# import Kratos
import KratosMultiphysics

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests to test_classes to create the suits
from chimera_analysis_test import FlowOverCrossFractionalStep
from chimera_analysis_test import FlowOverCrossMonolithic
from rotate_region_test import ChimeraRotateRegionTest

def AssembleTestSuites():
    ''' Populates the test suites to run.

    Populates the test suites to run. At least, it should populate the suites:
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
    smallSuite = suites['small']
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([ChimeraRotateRegionTest]))

    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)

    ## Validation suite. Big cases go here
    validationSuite = suites['validation']
    ### Single-Patch tests
    validationSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([FlowOverCrossMonolithic]))
    validationSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([FlowOverCrossFractionalStep]))
    ### Multi-Patch tests
    # Create a test suit that contains all the tests from every testCase
    # in the list:
    allSuite = suites['all']
    allSuite.addTests(nightSuite)
    allSuite.addTests(validationSuite)


    return suites

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning python tests ...")
    KratosUnittest.runTests(AssembleTestSuites())
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished python tests!")
