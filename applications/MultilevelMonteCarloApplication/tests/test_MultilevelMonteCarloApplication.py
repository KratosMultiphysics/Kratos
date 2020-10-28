# import Kratos
import KratosMultiphysics
import KratosMultiphysics.MultilevelMonteCarloApplication as KratosMLMC
import KratosMultiphysics.MeshingApplication as MeshingApplication

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
from test_multilevel_montecarlo import KratosMultilevelMonteCarloGeneralTests
from test_tools import TestTools
from test_xmcAlgorithm import TestXMCAlgorithm
from momentEstimatorTests import MomentEstimatorTest
from momentEstimatorTests import CombinedMomentEstimatorTest
from momentEstimatorTests import MultiMomentEstimatorTest
from momentEstimatorTests import MultiCombinedMomentEstimatorTest

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
    smallSuite.addTest(MomentEstimatorTest('test_update'))
    smallSuite.addTest(MomentEstimatorTest('test_value'))
    smallSuite.addTest(CombinedMomentEstimatorTest('test_update'))
    smallSuite.addTest(CombinedMomentEstimatorTest('test_value'))
    smallSuite.addTest(TestTools('test_normalInverseCDF'))
    smallSuite.addTest(TestTools('test_returnInput'))
    smallSuite.addTest(MultiMomentEstimatorTest('test_estimation_deterministic'))
    smallSuite.addTest(MultiMomentEstimatorTest('test_update'))
    smallSuite.addTest(MultiMomentEstimatorTest('test_construction_isParallel'))
    smallSuite.addTest(MultiCombinedMomentEstimatorTest('test_updateD0'))
    smallSuite.addTest(MultiCombinedMomentEstimatorTest('test_estimationD0'))

    # Create a test suit with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)
    nightSuite.addTest(MultiMomentEstimatorTest('test_estimation_random'))
    nightSuite.addTest(TestXMCAlgorithm('test_mc_Kratos'))
    nightSuite.addTest(KratosMultilevelMonteCarloGeneralTests('testMonteCarlo'))
    if(hasattr(MeshingApplication,"MmgProcess2D")):
        nightSuite.addTest(KratosMultilevelMonteCarloGeneralTests('testMultilevelMonteCarlo'))
        nightSuite.addTest(TestXMCAlgorithm('test_mlmc_Kratos'))
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
