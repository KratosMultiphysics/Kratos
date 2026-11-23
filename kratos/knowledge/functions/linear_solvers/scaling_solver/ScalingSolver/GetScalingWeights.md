---
type: C++ Method
title: GetScalingWeights
resource: linear_solvers/scaling_solver.h#L394-L411
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/linear_solvers/scaling_solver/ScalingSolver/GS2weights_row
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/linear_solvers/scaling_solver/ScalingSolver/Solve
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void GetScalingWeights( const SparseMatrixType& A, VectorType& aux)`

# Calls

- [GS2weights_row](../../../../functions/linear_solvers/scaling_solver/ScalingSolver/GS2weights_row.md)

# Called by

- [Solve](../../../../functions/linear_solvers/scaling_solver/ScalingSolver/Solve.md)