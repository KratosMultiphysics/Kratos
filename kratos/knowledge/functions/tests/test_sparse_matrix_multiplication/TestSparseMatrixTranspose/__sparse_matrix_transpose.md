---
type: Python Method
title: __sparse_matrix_transpose
resource: tests/test_sparse_matrix_multiplication.py#L51-L64
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/matrix_market_interface/ReadMatrixMarketMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/sparse_matrix_multiplication_utility/TransposeMatrix
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_sparse_matrix_multiplication/TestSparseMatrixTranspose/test_sparse_matrix_transpose_small
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sparse_matrix_multiplication/TestSparseMatrixTranspose/test_sparse_matrix_transpose_full
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __sparse_matrix_transpose(self, file_name = "auxiliar_files_for_python_unittest/sparse_matrix_files/A.mm"): # Read the matrices`

# Calls

- [ReadMatrixMarketMatrix](../../../../functions/sources/matrix_market_interface/ReadMatrixMarketMatrix.md)
- [TransposeMatrix](../../../../functions/utilities/sparse_matrix_multiplication_utility/TransposeMatrix.md)

# Called by

- [test_sparse_matrix_transpose_small](../../../../functions/tests/test_sparse_matrix_multiplication/TestSparseMatrixTranspose/test_sparse_matrix_transpose_small.md)
- [test_sparse_matrix_transpose_full](../../../../functions/tests/test_sparse_matrix_multiplication/TestSparseMatrixTranspose/test_sparse_matrix_transpose_full.md)