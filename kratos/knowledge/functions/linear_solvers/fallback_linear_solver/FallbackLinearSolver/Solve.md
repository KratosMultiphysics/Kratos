---
type: C++ Method
title: Solve
resource: linear_solvers/fallback_linear_solver.h#L267-L295
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/GetCurrentSolver
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/UpdateSolverIndex
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool Solve( SparseMatrixType& rA, VectorType& rX, VectorType& rB ) override`

# Calls

- [GetCurrentSolver](../../../../functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/GetCurrentSolver.md)
- [UpdateSolverIndex](../../../../functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/UpdateSolverIndex.md)