---
type: C++ Function
title: GetDofPosition
resource: includes/node.h#L1057-L1067
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/elements/edge_based_gradient_recovery_element/EdgeBasedGradientRecoveryElement<2>/EquationIdVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/edge_based_gradient_recovery_element/EdgeBasedGradientRecoveryElement<3>/EquationIdVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/embedded_nodal_variable_calculation_element_simplex/EmbeddedNodalVariableCalculationElementSimplex<double>/EquationIdVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/embedded_nodal_variable_calculation_element_simplex/EmbeddedNodalVariableCalculationElementSimplex<array_1d<double, 3>>/EquationIdVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise/ResidualBasedEliminationBuilderAndSolverComponentwise/Build
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise/ResidualBasedEliminationBuilderAndSolverComponentwise/ConstructGraph
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise/ParallelConstructGraph
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_bdf_displacement_scheme/ResidualBasedBDFDisplacementScheme/InitializeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_bdf_displacement_scheme/ResidualBasedBDFDisplacementScheme/Predict
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_bossak_displacement_scheme/ResidualBasedBossakDisplacementScheme/Predict
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_bossak_displacement_scheme/ResidualBasedBossakDisplacementScheme/InitializeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_pseudo_static_displacement_scheme/ResidualBasedPseudoStaticDisplacementScheme/Predict
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_node/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`inline unsigned int GetDofPosition(TVariableType const& rDofVariable) const`

# Calls

- [GetVariable](../../../functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable.md)

# Called by

- [EquationIdVector](../../../functions/elements/edge_based_gradient_recovery_element/EdgeBasedGradientRecoveryElement<2>/EquationIdVector.md)
- [EquationIdVector](../../../functions/elements/edge_based_gradient_recovery_element/EdgeBasedGradientRecoveryElement<3>/EquationIdVector.md)
- [EquationIdVector](../../../functions/elements/embedded_nodal_variable_calculation_element_simplex/EmbeddedNodalVariableCalculationElementSimplex<double>/EquationIdVector.md)
- [EquationIdVector](../../../functions/elements/embedded_nodal_variable_calculation_element_simplex/EmbeddedNodalVariableCalculationElementSimplex<array_1d<double, 3>>/EquationIdVector.md)
- [Build](../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise/ResidualBasedEliminationBuilderAndSolverComponentwise/Build.md)
- [ConstructGraph](../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise/ResidualBasedEliminationBuilderAndSolverComponentwise/ConstructGraph.md)
- [ParallelConstructGraph](../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise/ParallelConstructGraph.md)
- [InitializeSolutionStep](../../../functions/solving_strategies/schemes/residual_based_bdf_displacement_scheme/ResidualBasedBDFDisplacementScheme/InitializeSolutionStep.md)
- [Predict](../../../functions/solving_strategies/schemes/residual_based_bdf_displacement_scheme/ResidualBasedBDFDisplacementScheme/Predict.md)
- [Predict](../../../functions/solving_strategies/schemes/residual_based_bossak_displacement_scheme/ResidualBasedBossakDisplacementScheme/Predict.md)
- [InitializeSolutionStep](../../../functions/solving_strategies/schemes/residual_based_bossak_displacement_scheme/ResidualBasedBossakDisplacementScheme/InitializeSolutionStep.md)
- [Predict](../../../functions/solving_strategies/schemes/residual_based_pseudo_static_displacement_scheme/ResidualBasedPseudoStaticDisplacementScheme/Predict.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/sources/test_node/KRATOS_TEST_CASE_IN_SUITE.md)