---
type: C++ Function
title: MatrixAdd
resource: utilities/sparse_matrix_multiplication_utility.h#L364-L481
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/Matrix/fill
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/sparse_matrix_multiplication_utility/SortRows
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/sparse_matrix_multiplication_utility/CreateSolutionMatrix
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_sparse_matrix_multiplication/TestSparseMatrixSum/__sparse_matrix_sum
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void MatrixAdd( AMatrix& A, const BMatrix& B, const double Factor = 1.0 )`

# Calls

- [fill](../../../functions/includes/amatrix_interface/Matrix/fill.md)
- [SortRows](../../../functions/utilities/sparse_matrix_multiplication_utility/SortRows.md)
- [CreateSolutionMatrix](../../../functions/utilities/sparse_matrix_multiplication_utility/CreateSolutionMatrix.md)

# Called by

- [__sparse_matrix_sum](../../../functions/tests/test_sparse_matrix_multiplication/TestSparseMatrixSum/__sparse_matrix_sum.md)