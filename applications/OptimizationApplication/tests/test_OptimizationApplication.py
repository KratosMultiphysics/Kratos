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
from test_optimization_info import TestOptimizationInfo
from test_optimization_utils import TestOptimizationUtils
from test_mass_response_function import TestMassResponseFunctionBeams
from test_mass_response_function import TestMassResponseFunctionShells
from test_mass_response_function import TestMassResponseFunctionSolids
from test_linear_strain_energy_response_function import TestLinearStrainEnergyResponseFunction
from test_container_variable_data_holder import TestHistoricalContainerVariableDataHolder
from test_container_variable_data_holder import TestNodalContainerVariableDataHolder
from test_container_variable_data_holder import TestConditionContainerVariableDataHolder
from test_container_variable_data_holder import TestElementContainerVariableDataHolder
from test_container_variable_data_holder import TestConditionPropertiesContainerVariableDataHolder
from test_container_variable_data_holder import TestElementPropertiesContainerVariableDataHolder
from test_collective_variable_data_holder import TestCollectiveVariableDataHolder
from test_model_part_controllers import TestMdpaModelPartController
from test_container_variable_data_holder_utils import TestContainerVariableDataHolderUtils
from test_material_properties_control import TestMaterialPropertiesControl
from test_shape_control import TestShapeControl
from test_vertex_morphing_mappers import TestVertexMorphingNodalContainerVariableDataMapper
from test_vertex_morphing_mappers import TestVertexMorphingConditionContainerVariableDataMapper
from test_vertex_morphing_mappers import TestVertexMorphingElementContainerVariableDataMapper
from test_mass_optimization import TestMassOptimization

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
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestOptimizationInfo]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestOptimizationUtils]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestContainerVariableDataHolderUtils]))

    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestMassResponseFunctionBeams]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestMassResponseFunctionShells]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestMassResponseFunctionSolids]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestLinearStrainEnergyResponseFunction]))

    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestHistoricalContainerVariableDataHolder]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestNodalContainerVariableDataHolder]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestConditionContainerVariableDataHolder]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestElementContainerVariableDataHolder]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestConditionPropertiesContainerVariableDataHolder]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestElementPropertiesContainerVariableDataHolder]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCollectiveVariableDataHolder]))

    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestMdpaModelPartController]))

    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestMaterialPropertiesControl]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestShapeControl]))

    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestVertexMorphingNodalContainerVariableDataMapper]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestVertexMorphingConditionContainerVariableDataMapper]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestVertexMorphingElementContainerVariableDataMapper]))

    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestMassOptimization]))

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

    # Creating a test suit that contains all tests:
    allSuite = suites['all']
    allSuite.addTests(validationSuite)

    return suites

# ==============================================================================
# Main
# ==============================================================================
if __name__ == '__main__':
    km.Tester.SetVerbosity(km.Tester.Verbosity.PROGRESS)
    KratosUnittest.runTests(AssembleTestSuites())

# ==============================================================================
