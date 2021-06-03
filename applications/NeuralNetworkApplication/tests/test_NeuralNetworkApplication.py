# import Kratos

import KratosMultiphysics
import KratosMultiphysics.NeuralNetworkApplication

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

# Tests for preprocessing
from test_masking import TestMaskZeros
from test_normalizing import TestNormalizing

# Test for modeling
from test_modeling import TestInitialLayers
from test_modeling import TestInnerLayers
from test_modeling import TestBlockLayers
# from test_loading_model import TestLoadingModel

# Tests for training
# from test_training_parameters import TestOptimizer
# from test_training_parameters import TestLossFunction
# from test_training_parameters import TestMetrics
# from test_training_parameters import TestCallbacks
# from test_training import TestTraining

# Test for postprocessing
# from test_saving import TestSaving
# from test_testing import TestTesting

# Tests for output
# from test_plotters import TestHistoryPlotter
# from test_plotters import TestPredictionPlotter

# Tests for tuning 

# from test_tuning import TestLoadingHypermodel
# from test_tuning import TestBuildingHypermodel
# from test_tuning import TestTuning


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
    # Preprocessing
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestMaskZeros]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestNormalizing]))
    # Modeling
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestInitialLayers]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestInnerLayers]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestBlockLayers]))
    # smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestLoadingModel]))
    # smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestModeling]))
    # Training
    # smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestOptimizer]))
    # smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestMetrics]))
    # smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestLossFunction]))
    # smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCallbacks]))
    # Postprocessing
    # smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestSaving]))
    # Plotters
    # smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestHistoryPlotter]))
    # smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestPredictionPlotter]))
    # Tuning
    # smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestLoadingHypermodel]))
    # smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestBuildingHypermodel]))

    # Create a test suit with the selected tests

    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)
    # Data generation
    
    # Preprocessing
    
    # Modeling

    # Training
    # nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestTraining]))
    # Postprocessing
    # nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestTesting]))
    # Plotters
 
    # Tuning
    # nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestTuning]))

    validationSuite = suites['validation']
    validationSuite.addTest(nightSuite)
    # Create a test suit that contains all the tests from every testCase
    # in the list:
    allSuite = suites['all']
    allSuite.addTests(nightSuite)

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
