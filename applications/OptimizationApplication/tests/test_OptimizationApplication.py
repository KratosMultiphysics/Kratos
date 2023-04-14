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
import optimization_test_factory
import symmetry_utilities_tests.symmetry_tests
import test_execution_policies
import test_optimization_info
import test_optimization_utils
import test_mass_response_function
import test_linear_strain_energy_response_function
import test_model_part_utils
import test_model_part_controllers
import test_container_expression_utils
import test_container_expression
import test_collective_expressions
import test_buffered_dict
import test_standardized_responses

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
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_execution_policies.TestExecutionPolicies]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([symmetry_utilities_tests.symmetry_tests.SymmetryUtilitiesTest]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_buffered_dict.TestBufferedDict]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_optimization_info.TestOptimizationInfo]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_optimization_utils.TestOptimizationUtils]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_model_part_utils.TestModelPartUtils]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_container_expression_utils.TestContainerExpressionUtils]))

    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_mass_response_function.TestMassResponseFunctionBeams]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_mass_response_function.TestMassResponseFunctionShells]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_mass_response_function.TestMassResponseFunctionSolids]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_linear_strain_energy_response_function.TestLinearStrainEnergyResponseFunction]))

    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_standardized_responses.TestStandardizedObjective]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_standardized_responses.TestStandardizedConstraint]))

    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_container_expression.TestConditionPropertiesExpression]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_container_expression.TestElementPropertiesExpression]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_collective_expressions.TestCollectiveExpressions]))

    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_model_part_controllers.TestMdpaModelPartController]))

    # Adding nightly tests (tests that take < 10min)
    nightSuite = suites['nightly']

    # Adding small tests to nightly tests
    nightSuite.addTests(smallSuite)

    # Adding validation tests
    validationSuite = suites['validation']
    validationSuite.addTests(nightSuite)
    validationSuite.addTest(optimization_test_factory.top_opt_test('test_execution'))
    validationSuite.addTest(optimization_test_factory.mat_opt_test('test_execution'))
    validationSuite.addTest(optimization_test_factory.shell_shape_opt_test('test_execution'))
    validationSuite.addTest(optimization_test_factory.shell_thick_opt_test('test_execution'))

    # Creating a test suit that contains all tests:
    allSuite = suites['all']
    allSuite.addTests(validationSuite)

    return suites

# ==============================================================================
# Main
# ==============================================================================
if __name__ == '__main__':
    km.Tester.SetVerbosity(km.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    KratosUnittest.runTests(AssembleTestSuites())

# ==============================================================================
