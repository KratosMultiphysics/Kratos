---
type: C++ Method
title: InternalSystemSolveWithPhysics
resource: solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h#L460-L490
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/flags/IsNot
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/SystemSolveWithPhysics
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_lagrange_multiplier/ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier/SystemSolveWithPhysics
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void InternalSystemSolveWithPhysics( TSystemMatrixType& rA, TSystemVectorType& rDx, TSystemVectorType& rb, ModelPart& rModelPart )`

# Calls

- [IsNot](../../../../../functions/containers/flags/IsNot.md)

# Called by

- [SystemSolveWithPhysics](../../../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/SystemSolveWithPhysics.md)
- [SystemSolveWithPhysics](../../../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_lagrange_multiplier/ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier/SystemSolveWithPhysics.md)