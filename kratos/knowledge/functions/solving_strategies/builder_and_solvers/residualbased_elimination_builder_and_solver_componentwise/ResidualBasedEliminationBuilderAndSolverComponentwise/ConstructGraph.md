---
type: C++ Method
title: ConstructGraph
resource: solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise.h#L543-L602
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/node/GetDofPosition
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise/ResidualBasedEliminationBuilderAndSolverComponentwise/ResizeAndInitializeVectors
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ConstructGraph(TSystemMatrixType& A)`

# Calls

- [GetDofPosition](../../../../../functions/includes/node/GetDofPosition.md)

# Called by

- [ResizeAndInitializeVectors](../../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise/ResidualBasedEliminationBuilderAndSolverComponentwise/ResizeAndInitializeVectors.md)