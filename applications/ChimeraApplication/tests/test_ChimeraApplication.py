# import Kratos
import KratosMultiphysics

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests to test_classes to create the suits
#sys.path.insert(0, 'chimera_monolithic_simple_test')
#from chimera_analysis_test import FlowOverCylinderMonolithic
#from chimera_analysis_test import FlowOverCylinderFractionalStep
from chimera_analysis_test import FlowOverCrossFractionalStep
from chimera_analysis_test import FlowOverCrossMonolithic
from rotate_region_test import ChimeraRotateRegionTest
#from chimera_analysis_test import MonolithicMultiPatch
#from chimera_analysis_test import FractionalStepMultiPatch

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
    # - testSmallExample
    smallSuite = suites['small']
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([ChimeraRotateRegionTest]))
    ### Single-Patch tests
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([FlowOverCrossMonolithic]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([FlowOverCrossFractionalStep]))
    ### Multi-Patch tests
    #smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([MonolithicMultiPatch]))
    #smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([FractionalStepMultiPatch]))

    # Create a test suit with the selected tests
    # nightSuite will contain the following tests:
    # - testSmallExample
    # - testNightlyFirstExample
    # - testNightlySecondExample
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)

    ## Validation suite. Big cases go here
    #validationSuite = suites['validation']
    #validationSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([FlowOverCylinderMonolithic]))
    #validationSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([FlowOverCylinderFractionalStep]))
    # Create a test suit that contains all the tests from every testCase
    # in the list:
    allSuite = suites['all']
    allSuite.addTests(nightSuite)


    return suites

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning python tests ...")
    KratosUnittest.runTests(AssembleTestSuites())
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished python tests!")
