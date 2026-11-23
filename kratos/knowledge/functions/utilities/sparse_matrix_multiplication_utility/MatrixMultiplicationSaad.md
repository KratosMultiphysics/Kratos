---
type: C++ Function
title: MatrixMultiplicationSaad
resource: utilities/sparse_matrix_multiplication_utility.h#L133-L237
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/sparse_matrix_multiplication_utility/SortRows
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/sparse_matrix_multiplication_utility/CreateSolutionMatrix
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_sparse_matrix_multiplication/TestSparseMatrixMultiplication/__sparse_matrix_multiplication
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/sparse_matrix_multiplication_utility/MatrixMultiplication
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void MatrixMultiplicationSaad( const AMatrix& A, const BMatrix& B, CMatrix& C )`

# Calls

- [SortRows](../../../functions/utilities/sparse_matrix_multiplication_utility/SortRows.md)
- [CreateSolutionMatrix](../../../functions/utilities/sparse_matrix_multiplication_utility/CreateSolutionMatrix.md)

# Called by

- [__sparse_matrix_multiplication](../../../functions/tests/test_sparse_matrix_multiplication/TestSparseMatrixMultiplication/__sparse_matrix_multiplication.md)
- [MatrixMultiplication](../../../functions/utilities/sparse_matrix_multiplication_utility/MatrixMultiplication.md)