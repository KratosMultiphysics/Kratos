# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
from generalTests import KratosNeuralNetworkGeneralTests
# Tests for data generation
from test_data_generation import TestDataGeneration1I1OVariablesText
from test_data_generation import TestDataGeneration1I2OVariablesText
from test_data_generation import TestDataGeneration1IVectorOVariablesText
# from test_data_generation import TestDataGeneration1I1OVariablesHDF5
# from test_data_generation import TestDataGeneration1I2OVariablesHDF5
# from test_data_generation import TestDataGeneration1IVectorOVariablesHDF5
# from test_data_generation import TestDataGenerationGlobals
# from test_data_generation import TestDataGenerationGeometrical

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
    smallSuite.addTest(KratosNeuralNetworkGeneralTests('testSmallExample'))
    # Data generation
    # nightSuite.addTest(TestDataGeneration1I1OVariablesHDF5('test_execution'))
    smallSuite.addTest(TestDataGeneration1I1OVariablesText('test_execution'))
    # nightSuite.addTest(TestDataGeneration1I2OVariablesHDF5('test_execution'))
    smallSuite.addTest(TestDataGeneration1I2OVariablesText('test_execution'))
    # nightSuite.addTest(TestDataGeneration1IVectorOVariablesHDF5('test_execution'))
    smallSuite.addTest(TestDataGeneration1IVectorOVariablesText('test_execution'))

    # Create a test suit with the selected tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)
    # Data generation

    validationSuite = suites['validation']
    validationSuite.addTest(nightSuite)
    # Create a test suit that contains all the tests from every testCase
    # in the list:
    allSuite = suites['all']
    allSuite.addTests(nightSuite)

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
