---
type: C++ Method
title: ProvideAdditionalData
resource: linear_solvers/fallback_linear_solver.h#L378-L393
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/GetCurrentSolver
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ProvideAdditionalData( SparseMatrixType& rA, VectorType& rX, VectorType& rB, typename ModelPart::DofsArrayType& rDoFSet, ModelPart& rModelPart ) override`

# Calls

- [GetCurrentSolver](../../../../functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/GetCurrentSolver.md)