# ==============================================================================
# Imports
# ==============================================================================


# Import Kratos "wrapper" for unittests
import KratosMultiphysics as km
import KratosMultiphysics.KratosUnittest as KratosUnittest

# ==============================================================================
# Import the tests or test_classes to create the suits
# ==============================================================================

# Small tests
from optimization_test_factory import top_opt_test
from optimization_test_factory import mat_opt_test
from optimization_test_factory import shell_shape_opt_test
from optimization_test_factory import shell_thick_opt_test
from test_execution_policies import TestExecutionPolicies
from test_optimization_variable_utils import TestOptimizationVariableUtils
from test_mass_response_function import TestMassResponseFunctionBeams
from test_mass_response_function import TestMassResponseFunctionShells
from test_mass_response_function import TestMassResponseFunctionSolids
from test_entity_specific_properties_process import TestEntitySpecificPropertiesProcess
from test_linear_strain_energy_response_function import TestLinearStrainEnergyResponseFunctionBase
from test_material_properties_control import TestMaterialPropertiesControl
from test_mass_optimization import TestMassOptimization
from test_shape_control import TestShapeControl
from test_optimization_info import TestOptimizationInfo
from test_container_data import TestHistoricalContainerVariableDataHolder
from test_container_data import TestNodalContainerVariableDataHolder
from test_container_data import TestConditionContainerVariableDataHolder
from test_container_data import TestElementContainerVariableDataHolder
from test_container_data import TestConditionPropertiesContainerVariableDataHolder
from test_container_data import TestElementPropertiesContainerVariableDataHolder
from test_container_data_utils import TestContainerVariableDataHolderUtils

# Nightly tests

# Validation tests

# ==============================================================================
# Test assembly
# ==============================================================================
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

    # Adding small tests (tests that take < 1s)
    smallSuite = suites['small']
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestExecutionPolicies]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestOptimizationVariableUtils]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestMassResponseFunctionBeams]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestMassResponseFunctionShells]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestMassResponseFunctionSolids]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestEntitySpecificPropertiesProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestLinearStrainEnergyResponseFunctionBase]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestMaterialPropertiesControl]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestShapeControl]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestOptimizationInfo]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestHistoricalContainerVariableDataHolder]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestNodalContainerVariableDataHolder]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestConditionContainerVariableDataHolder]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestElementContainerVariableDataHolder]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestConditionPropertiesContainerVariableDataHolder]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestElementPropertiesContainerVariableDataHolder]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestContainerVariableDataHolderUtils]))

    # Adding nightly tests (tests that take < 10min)
    nightSuite = suites['nightly']

    # Adding small tests to nightly tests
    nightSuite.addTests(smallSuite)

    # Adding validation tests
    validationSuite = suites['validation']
    validationSuite.addTests(nightSuite)
    validationSuite.addTest(top_opt_test('test_execution'))
    validationSuite.addTest(mat_opt_test('test_execution'))
    validationSuite.addTest(shell_shape_opt_test('test_execution'))
    validationSuite.addTest(shell_thick_opt_test('test_execution'))
    validationSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestMassOptimization]))

    # Creating a test suit that contains all tests:
    allSuite = suites['all']
    allSuite.addTests(validationSuite)

    return suites

# ==============================================================================
# Main
# ==============================================================================
if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())

# ==============================================================================
