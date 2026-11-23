---
type: C++ Method
title: SymmetricScaling
resource: linear_solvers/scaling_solver.h#L352-L371
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/linear_solvers/scaling_solver/ScalingSolver/perform_matrix_scaling_row
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/linear_solvers/scaling_solver/ScalingSolver/Solve
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void SymmetricScaling( SparseMatrixType& A, const VectorType& aux)`

# Calls

- [perform_matrix_scaling_row](../../../../functions/linear_solvers/scaling_solver/ScalingSolver/perform_matrix_scaling_row.md)

# Called by

- [Solve](../../../../functions/linear_solvers/scaling_solver/ScalingSolver/Solve.md)