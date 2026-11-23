---
type: C++ Function
title: MatrixMultiplicationRMerge
resource: utilities/sparse_matrix_multiplication_utility.h#L247-L356
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/sparse_matrix_multiplication_utility/ProdRowWidth
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/sparse_matrix_multiplication_utility/ProdRow
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

`static void MatrixMultiplicationRMerge( const AMatrix &A, const BMatrix &B, CMatrix &C )`

# Calls

- [ProdRowWidth](../../../functions/utilities/sparse_matrix_multiplication_utility/ProdRowWidth.md)
- [ProdRow](../../../functions/utilities/sparse_matrix_multiplication_utility/ProdRow.md)
- [CreateSolutionMatrix](../../../functions/utilities/sparse_matrix_multiplication_utility/CreateSolutionMatrix.md)

# Called by

- [__sparse_matrix_multiplication](../../../functions/tests/test_sparse_matrix_multiplication/TestSparseMatrixMultiplication/__sparse_matrix_multiplication.md)
- [MatrixMultiplication](../../../functions/utilities/sparse_matrix_multiplication_utility/MatrixMultiplication.md)