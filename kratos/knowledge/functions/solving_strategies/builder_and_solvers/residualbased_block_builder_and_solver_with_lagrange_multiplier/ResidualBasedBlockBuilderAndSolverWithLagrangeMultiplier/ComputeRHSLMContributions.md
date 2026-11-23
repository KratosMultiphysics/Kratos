---
type: C++ Method
title: ComputeRHSLMContributions
resource: solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_lagrange_multiplier.h#L1116-L1168
visibility: private
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
---

# Signature

`void ComputeRHSLMContributions( TSystemVectorType& rbLM, const double ScaleFactor = 1.0 )`

# Called by

- [ApplyRHSConstraints](../../../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_lagrange_multiplier/ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier/ApplyRHSConstraints.md)
- [ApplyConstraints](../../../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_lagrange_multiplier/ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier/ApplyConstraints.md)