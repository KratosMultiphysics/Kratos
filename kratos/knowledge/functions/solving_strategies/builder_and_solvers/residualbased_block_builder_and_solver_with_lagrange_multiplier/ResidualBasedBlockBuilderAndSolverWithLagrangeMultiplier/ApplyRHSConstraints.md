---
type: C++ Method
title: ApplyRHSConstraints
resource: solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_lagrange_multiplier.h#L492-L555
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/spaces/ublas_space/GetMaxDiagonal
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spaces/ublas_space/GetDiagonalNorm
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/flags/IsNot
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_lagrange_multiplier/ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier/ComputeRHSLMContributions
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ApplyRHSConstraints( typename TSchemeType::Pointer pScheme, ModelPart& rModelPart, TSystemVectorType& rb ) override`

# Calls

- [GetMaxDiagonal](../../../../../functions/spaces/ublas_space/GetMaxDiagonal.md)
- [GetDiagonalNorm](../../../../../functions/spaces/ublas_space/GetDiagonalNorm.md)
- [IsNot](../../../../../functions/containers/flags/IsNot.md)
- [ComputeRHSLMContributions](../../../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_lagrange_multiplier/ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier/ComputeRHSLMContributions.md)