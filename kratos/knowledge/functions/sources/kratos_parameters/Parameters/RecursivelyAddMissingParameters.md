---
type: C++ Method
title: RecursivelyAddMissingParameters
resource: sources/kratos_parameters.cpp#L1496-L1520
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/IsSubParameter
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/pointer_vector_map_iterator/PointerVectorMapIterator/key
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/future/solving_strategies/schemes/static_scheme/StaticScheme/GetDefaultParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/strategies/implicit_strategy/ImplicitStrategy/GetDefaultParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/strategies/linear_strategy/LinearStrategy/GetDefaultParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/apply_ray_casting_interface_recognition_process/ApplyRayCastingInterfaceRecognitionProcess<TDim>/GetDefaultParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_distance_to_skin_process/CalculateDistanceToSkinProcess<TDim>/GetDefaultParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/modelers/import_obj_modeler/ImportOBJModeler/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/GetDefaultParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/GetDefaultParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_lagrange_multiplier/ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier/GetDefaultParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver/ResidualBasedEliminationBuilderAndSolver/GetDefaultParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/GetDefaultParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/convergencecriterias/and_criteria/And_Criteria/GetDefaultParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/convergencecriterias/displacement_criteria/DisplacementCriteria/GetDefaultParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/convergencecriterias/or_criteria/Or_Criteria/GetDefaultParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/convergencecriterias/residual_criteria/ResidualCriteria/GetDefaultParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_bdf_custom_scheme/ResidualBasedBDFCustomScheme/GetDefaultParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_bdf_displacement_scheme/ResidualBasedBDFDisplacementScheme/GetDefaultParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_bdf_scheme/ResidualBasedBDFScheme/GetDefaultParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_bossak_displacement_scheme/ResidualBasedBossakDisplacementScheme/GetDefaultParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_implicit_time_scheme/ResidualBasedImplicitTimeScheme/GetDefaultParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_newmark_displacement_scheme/ResidualBasedNewmarkDisplacementScheme/GetDefaultParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_pseudo_static_displacement_scheme/ResidualBasedPseudoStaticDisplacementScheme/GetDefaultParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residualbased_incrementalupdate_static_scheme/ResidualBasedIncrementalUpdateStaticScheme/GetDefaultParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/adaptive_residualbased_newton_raphson_strategy/AdaptiveResidualBasedNewtonRaphsonStrategy/GetDefaultParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/arc_length_strategy/ArcLengthStrategy/GetDefaultParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/explicit_solving_strategy/ExplicitSolvingStrategy/GetDefaultParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/explicit_solving_strategy_bfecc/ExplicitSolvingStrategyBFECC/GetDefaultParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/explicit_solving_strategy_runge_kutta/ExplicitSolvingStrategyRungeKutta/GetDefaultParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/implicit_solving_strategy/ImplicitSolvingStrategy/GetDefaultParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/line_search_strategy/GetDefaultParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/residualbased_linear_strategy/ResidualBasedLinearStrategy/GetDefaultParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/residualbased_newton_raphson_strategy/GetDefaultParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/RecursivelyValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_kratos_parameters/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_recursively_add_missing_parameters
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Parameters::RecursivelyAddMissingParameters(const Parameters& rDefaultParameters)`

# Calls

- [IsSubParameter](../../../../functions/sources/kratos_parameters/Parameters/IsSubParameter.md)
- [key](../../../../functions/containers/pointer_vector_map_iterator/PointerVectorMapIterator/key.md)

# Called by

- [GetDefaultParameters](../../../../functions/future/solving_strategies/schemes/static_scheme/StaticScheme/GetDefaultParameters.md)
- [GetDefaultParameters](../../../../functions/future/solving_strategies/strategies/implicit_strategy/ImplicitStrategy/GetDefaultParameters.md)
- [GetDefaultParameters](../../../../functions/future/solving_strategies/strategies/linear_strategy/LinearStrategy/GetDefaultParameters.md)
- [GetDefaultParameters](../../../../functions/processes/apply_ray_casting_interface_recognition_process/ApplyRayCastingInterfaceRecognitionProcess<TDim>/GetDefaultParameters.md)
- [GetDefaultParameters](../../../../functions/processes/calculate_distance_to_skin_process/CalculateDistanceToSkinProcess<TDim>/GetDefaultParameters.md)
- [__init__](../../../../functions/python_scripts/modelers/import_obj_modeler/ImportOBJModeler/__init__.md)
- [GetDefaultParameters](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/GetDefaultParameters.md)
- [GetDefaultParameters](../../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/GetDefaultParameters.md)
- [GetDefaultParameters](../../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_lagrange_multiplier/ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier/GetDefaultParameters.md)
- [GetDefaultParameters](../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver/ResidualBasedEliminationBuilderAndSolver/GetDefaultParameters.md)
- [GetDefaultParameters](../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/GetDefaultParameters.md)
- [GetDefaultParameters](../../../../functions/solving_strategies/convergencecriterias/and_criteria/And_Criteria/GetDefaultParameters.md)
- [GetDefaultParameters](../../../../functions/solving_strategies/convergencecriterias/displacement_criteria/DisplacementCriteria/GetDefaultParameters.md)
- [GetDefaultParameters](../../../../functions/solving_strategies/convergencecriterias/or_criteria/Or_Criteria/GetDefaultParameters.md)
- [GetDefaultParameters](../../../../functions/solving_strategies/convergencecriterias/residual_criteria/ResidualCriteria/GetDefaultParameters.md)
- [GetDefaultParameters](../../../../functions/solving_strategies/schemes/residual_based_bdf_custom_scheme/ResidualBasedBDFCustomScheme/GetDefaultParameters.md)
- [GetDefaultParameters](../../../../functions/solving_strategies/schemes/residual_based_bdf_displacement_scheme/ResidualBasedBDFDisplacementScheme/GetDefaultParameters.md)
- [GetDefaultParameters](../../../../functions/solving_strategies/schemes/residual_based_bdf_scheme/ResidualBasedBDFScheme/GetDefaultParameters.md)
- [GetDefaultParameters](../../../../functions/solving_strategies/schemes/residual_based_bossak_displacement_scheme/ResidualBasedBossakDisplacementScheme/GetDefaultParameters.md)
- [GetDefaultParameters](../../../../functions/solving_strategies/schemes/residual_based_implicit_time_scheme/ResidualBasedImplicitTimeScheme/GetDefaultParameters.md)
- [GetDefaultParameters](../../../../functions/solving_strategies/schemes/residual_based_newmark_displacement_scheme/ResidualBasedNewmarkDisplacementScheme/GetDefaultParameters.md)
- [GetDefaultParameters](../../../../functions/solving_strategies/schemes/residual_based_pseudo_static_displacement_scheme/ResidualBasedPseudoStaticDisplacementScheme/GetDefaultParameters.md)
- [GetDefaultParameters](../../../../functions/solving_strategies/schemes/residualbased_incrementalupdate_static_scheme/ResidualBasedIncrementalUpdateStaticScheme/GetDefaultParameters.md)
- [GetDefaultParameters](../../../../functions/solving_strategies/strategies/adaptive_residualbased_newton_raphson_strategy/AdaptiveResidualBasedNewtonRaphsonStrategy/GetDefaultParameters.md)
- [GetDefaultParameters](../../../../functions/solving_strategies/strategies/arc_length_strategy/ArcLengthStrategy/GetDefaultParameters.md)
- [GetDefaultParameters](../../../../functions/solving_strategies/strategies/explicit_solving_strategy/ExplicitSolvingStrategy/GetDefaultParameters.md)
- [GetDefaultParameters](../../../../functions/solving_strategies/strategies/explicit_solving_strategy_bfecc/ExplicitSolvingStrategyBFECC/GetDefaultParameters.md)
- [GetDefaultParameters](../../../../functions/solving_strategies/strategies/explicit_solving_strategy_runge_kutta/ExplicitSolvingStrategyRungeKutta/GetDefaultParameters.md)
- [GetDefaultParameters](../../../../functions/solving_strategies/strategies/implicit_solving_strategy/ImplicitSolvingStrategy/GetDefaultParameters.md)
- [GetDefaultParameters](../../../../functions/solving_strategies/strategies/line_search_strategy/GetDefaultParameters.md)
- [GetDefaultParameters](../../../../functions/solving_strategies/strategies/residualbased_linear_strategy/ResidualBasedLinearStrategy/GetDefaultParameters.md)
- [GetDefaultParameters](../../../../functions/solving_strategies/strategies/residualbased_newton_raphson_strategy/GetDefaultParameters.md)
- [RecursivelyValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/RecursivelyValidateAndAssignDefaults.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_kratos_parameters/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_recursively_add_missing_parameters](../../../../functions/tests/test_kratos_parameters/TestParameters/test_recursively_add_missing_parameters.md)