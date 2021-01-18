# import Kratos
import KratosMultiphysics
import KratosMultiphysics.MORApplication

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
from generalTests import KratosMORGeneralTests
from test_mor_utilities import TestMORUtilities
from test_irka_strategies import TestIrkaStrategies
from test_acoustic_elements import TestAcousticElements
from test_acoustic_conditions import TestAcousticConditions

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
    smallSuite = suites['small']
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestMORUtilities]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestIrkaStrategies]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAcousticElements]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAcousticConditions]))

    # Create a test suit with the selected tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)

    # Create a test suit that contains all the tests from every testCase
    allSuite = suites['all']
    allSuite.addTests(nightSuite)

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
