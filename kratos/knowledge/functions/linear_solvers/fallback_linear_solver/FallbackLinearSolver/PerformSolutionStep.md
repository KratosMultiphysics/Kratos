---
type: C++ Method
title: PerformSolutionStep
resource: linear_solvers/fallback_linear_solver.h#L231-L238
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/GetCurrentSolver
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool PerformSolutionStep( SparseMatrixType& rA, VectorType& rX, VectorType& rB ) override`

# Calls

- [GetCurrentSolver](../../../../functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/GetCurrentSolver.md)