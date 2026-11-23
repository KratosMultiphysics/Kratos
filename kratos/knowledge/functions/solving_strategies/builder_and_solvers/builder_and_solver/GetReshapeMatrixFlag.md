---
type: C++ Function
title: GetReshapeMatrixFlag
resource: solving_strategies/builder_and_solvers/builder_and_solver.h#L220-L223
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/ResizeAndInitializeVectors
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/ResizeAndInitializeVectors
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver/ResidualBasedEliminationBuilderAndSolver/ResizeAndInitializeVectors
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise/ResidualBasedEliminationBuilderAndSolverComponentwise/ResizeAndInitializeVectors
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/ResizeAndInitializeVectors
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool GetReshapeMatrixFlag() const`

# Called by

- [ResizeAndInitializeVectors](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/ResizeAndInitializeVectors.md)
- [ResizeAndInitializeVectors](../../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/ResizeAndInitializeVectors.md)
- [ResizeAndInitializeVectors](../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver/ResidualBasedEliminationBuilderAndSolver/ResizeAndInitializeVectors.md)
- [ResizeAndInitializeVectors](../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise/ResidualBasedEliminationBuilderAndSolverComponentwise/ResizeAndInitializeVectors.md)
- [ResizeAndInitializeVectors](../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/ResizeAndInitializeVectors.md)