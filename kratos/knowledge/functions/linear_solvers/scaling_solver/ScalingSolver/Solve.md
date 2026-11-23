---
type: C++ Method
title: Solve
resource: linear_solvers/scaling_solver.h#L224-L265
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/linear_solvers/scaling_solver/ScalingSolver/GetScalingWeights
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/linear_solvers/scaling_solver/ScalingSolver/SymmetricScaling
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override`

# Calls

- [GetScalingWeights](../../../../functions/linear_solvers/scaling_solver/ScalingSolver/GetScalingWeights.md)
- [SymmetricScaling](../../../../functions/linear_solvers/scaling_solver/ScalingSolver/SymmetricScaling.md)