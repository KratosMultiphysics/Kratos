---
type: C++ Function
title: GetDiagonalNorm
resource: spaces/ublas_space.h#L807-L825
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_lagrange_multiplier/ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier/ApplyRHSConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_lagrange_multiplier/ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier/ApplyConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spaces/ublas_space/GetScaleNorm
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static double GetDiagonalNorm(const MatrixType& rA)`

# Called by

- [ApplyRHSConstraints](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_lagrange_multiplier/ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier/ApplyRHSConstraints.md)
- [ApplyConstraints](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_lagrange_multiplier/ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier/ApplyConstraints.md)
- [GetScaleNorm](../../../functions/spaces/ublas_space/GetScaleNorm.md)