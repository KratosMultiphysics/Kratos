---
type: C++ Method
title: IsNotConsistent
resource: linear_solvers/fallback_linear_solver.h#L497-L504
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/GetCurrentSolver
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool IsNotConsistent( SparseMatrixType& rA, VectorType& rX, VectorType& rB ) override`

# Calls

- [GetCurrentSolver](../../../../functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/GetCurrentSolver.md)