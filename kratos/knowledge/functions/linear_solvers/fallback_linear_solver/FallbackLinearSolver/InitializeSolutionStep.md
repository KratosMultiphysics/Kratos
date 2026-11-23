---
type: C++ Method
title: InitializeSolutionStep
resource: linear_solvers/fallback_linear_solver.h#L205-L212
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/GetCurrentSolver
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void InitializeSolutionStep( SparseMatrixType& rA, VectorType& rX, VectorType& rB ) override`

# Calls

- [GetCurrentSolver](../../../../functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/GetCurrentSolver.md)