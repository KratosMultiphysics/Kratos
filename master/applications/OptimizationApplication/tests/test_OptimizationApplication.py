# ==============================================================================
# Imports
# ==============================================================================


# Import Kratos "wrapper" for unittests
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

# ==============================================================================
# Import the tests or test_classes to create the suits
# ==============================================================================

# Small tests
import algorithm_tests.analysis_based_tests.algorithm_steepest_descent_qnbb
import algorithm_tests.analysis_based_tests.algorithm_steepest_descent_qnbb.test_steepest_descent_analysis
import optimization_test_factory
import symmetry_utilities_tests.symmetry_tests
import test_execution_policies
import test_optimization_info
import test_optimization_utils
import test_response_utilities
import responses_tests.test_response_routine
import responses_tests.test_overhang_response_function
import responses_tests.test_mass_response_function
import responses_tests.test_linear_strain_energy_response_function
import responses_tests.test_standardized_responses
import responses_tests.test_geometric_centroid_deviation_response_function
import responses_tests.test_combined_response_function
import responses_tests.test_discrete_value_residual_response_function
import test_model_part_utils
import test_model_part_controllers
import test_connectivity_preserving_model_part_controller
import test_container_expression_utils
import test_container_expression
import test_collective_expressions
import test_sigmoidal_projection
import test_buffered_dict
import control.test_master_control
import control.material.test_material_properties_control
import control.thickness.test_shell_thickness_control
import control.shape.test_vm_shape_control
import control.material.test_simp_control
import filtering.implicit_filters_tests
import filtering.explicit_filters_tests
import test_component_data_view
import process_tests.test_optimization_problem_vtu_output_process
import process_tests.test_optimization_problem_ascii_output_process
import process_tests.test_optimization_problem_graph_output_process
import algorithm_tests.test_convergence
import algorithm_tests.test_line_search
import algorithm_tests.test_algorithm_steepest_descent
import algorithm_tests.analysis_based_tests.algorithm_steepest_descent_qnbb.test_steepest_descent_analysis
import algorithm_tests.analysis_based_tests.algorithm_steepest_descent.test_steepest_descent_analysis
import algorithm_tests.analysis_based_tests.algorithm_nesterov_accelerated_gradient.test_nestervo_accelerated_gradient_analysis
import algorithm_tests.analysis_based_tests.algorithm_gradient_projection.test_gradient_projection
import algorithm_tests.analysis_based_tests.algorithm_relaxed_gradient_projection.test_relaxed_gradient_projection
import algorithm_tests.nlopt_tests.mma_shell_thickness_opt.test_mma_optimizer

# Nightly tests

# Validation tests

# ==============================================================================
# Test assembly
# ==============================================================================
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

    # Adding small tests (tests that take < 1s)
    smallSuite = suites['small']
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_execution_policies.TestExecutionPolicies]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([symmetry_utilities_tests.symmetry_tests.SymmetryUtilitiesTest]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_buffered_dict.TestBufferedDict]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_component_data_view.TestComponentDataView]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_optimization_info.TestOptimizationInfo]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_optimization_utils.TestOptimizationUtils]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_model_part_utils.TestOptAppModelPartUtils]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_model_part_utils.TestModelPartUtilities]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_response_utilities.TestResponseUtilities]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_container_expression_utils.TestContainerExpressionUtils]))

    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([responses_tests.test_response_routine.TestResponseRoutine]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([responses_tests.test_standardized_responses.TestStandardizedObjective]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([responses_tests.test_standardized_responses.TestStandardizedConstraint]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([responses_tests.test_mass_response_function.TestMassResponseFunctionBeams]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([responses_tests.test_mass_response_function.TestMassResponseFunctionShells]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([responses_tests.test_mass_response_function.TestMassResponseFunctionSolids]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([responses_tests.test_linear_strain_energy_response_function.TestLinearStrainEnergyResponseFunction]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([responses_tests.test_overhang_response_function.TestOverHangResponseFunction]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([responses_tests.test_geometric_centroid_deviation_response_function.TestGeometricCentroidDeviationResponseFunction]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([responses_tests.test_combined_response_function.TestCombinedResponseFunction]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([responses_tests.test_discrete_value_residual_response_function.TestDiscreteValueResidualResponseFunctionExact]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([responses_tests.test_discrete_value_residual_response_function.TestDiscreteValueResidualResponseFunctionLogarithm]))

    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_container_expression.TestConditionPropertiesExpression]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_container_expression.TestElementPropertiesExpression]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_collective_expressions.TestCollectiveExpressions]))

    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_sigmoidal_projection.TestSigmoidalProjection]))

    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_model_part_controllers.TestMdpaModelPartController]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_connectivity_preserving_model_part_controller.TestConnectivityPreservingModelPartController]))

    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([control.test_master_control.TestMassterControl]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([control.material.test_material_properties_control.TestMaterialPropertiesControl]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([control.material.test_simp_control.TestSimpControl]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([control.thickness.test_shell_thickness_control.TestShellThicknessControl]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([control.shape.test_vm_shape_control.TestVMShapeControlShell]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([control.shape.test_vm_shape_control.TestVMShapeControlSolid]))

    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([filtering.implicit_filters_tests.HelmholtzAnalysisTest]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([filtering.explicit_filters_tests.TestExplicitFilterConsistency]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([filtering.explicit_filters_tests.TestExplicitFilterReference]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([process_tests.test_optimization_problem_vtu_output_process.TestOptimizationProblemVtuOutputProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([process_tests.test_optimization_problem_ascii_output_process.TestOptimizationProblemAsciiOutputProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([process_tests.test_optimization_problem_graph_output_process.TestOptimizationProblemGraphOutputProcess]))

    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([algorithm_tests.test_convergence.TestConvergence]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([algorithm_tests.test_line_search.TestLineSearch]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([algorithm_tests.test_algorithm_steepest_descent.TestAlgorithmSteepestDescent]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([algorithm_tests.analysis_based_tests.algorithm_steepest_descent.test_steepest_descent_analysis.TestSteepestDescentAnalysis]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([algorithm_tests.analysis_based_tests.algorithm_steepest_descent_qnbb.test_steepest_descent_analysis.TestQNBBSteepestDescentAnalysis]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([algorithm_tests.nlopt_tests.mma_shell_thickness_opt.test_mma_optimizer.TestNLOPTOptimizer]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([algorithm_tests.analysis_based_tests.algorithm_gradient_projection.test_gradient_projection.TestGradientProjectionAnalysis]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([algorithm_tests.analysis_based_tests.algorithm_relaxed_gradient_projection.test_relaxed_gradient_projection.TestRelaxedGradientProjectionAnalysis]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([algorithm_tests.analysis_based_tests.algorithm_nesterov_accelerated_gradient.test_nestervo_accelerated_gradient_analysis.TestNesterovAcceleratedGradientAnalysis]))



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
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.runTests(AssembleTestSuites())

# ==============================================================================
