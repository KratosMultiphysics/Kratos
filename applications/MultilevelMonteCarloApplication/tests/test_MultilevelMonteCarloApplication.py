# import Kratos
import KratosMultiphysics
import KratosMultiphysics.MultilevelMonteCarloApplication as KratosMLMC
import KratosMultiphysics.MeshingApplication as MeshingApplication

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
from test_multilevel_montecarlo import KratosMultilevelMonteCarloGeneralTests
from test_momentEstimator import TestMomentEstimator
from test_momentEstimator import TestCombinedMomentEstimator
from test_tools import TestTools
from test_xmcAlgorithm import TestXMCAlgorithm

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

    # Create a test suit with the selected tests
    smallSuite = suites['small']
    smallSuite.addTest(TestMomentEstimator('test_update'))
    smallSuite.addTest(TestMomentEstimator('test_value'))
    smallSuite.addTest(TestCombinedMomentEstimator('test_update'))
    smallSuite.addTest(TestCombinedMomentEstimator('test_value'))
    smallSuite.addTest(TestTools('test_normalInverseCDF'))

    # Create a test suit with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)
    nightSuite.addTest(TestXMCAlgorithm('test_mc_asynchronous_Kratos'))
    nightSuite.addTest(KratosMultilevelMonteCarloGeneralTests('testMonteCarlo'))
    if(hasattr(MeshingApplication,"MmgProcess2D")):
        nightSuite.addTest(KratosMultilevelMonteCarloGeneralTests('testMultilevelMonteCarlo'))
        nightSuite.addTest(TestXMCAlgorithm('test_mlmc_asynchronous_Kratos'))
    else:
        print("MMG process is not compiled and the corresponding tests will not be executed")

    # Create a test suit that contains all the tests
    allSuite = suites['all']
    allSuite.addTests(nightSuite)
    
    return suites

if __name__ == '__main__':
    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning python tests ...")
    KratosUnittest.runTests(AssembleTestSuites())
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished python tests!")
