---
type: Python Method
title: __sparse_matrix_sum
resource: tests/test_sparse_matrix_multiplication.py#L21-L39
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/matrix_market_interface/ReadMatrixMarketMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/sparse_matrix_multiplication_utility/MatrixAdd
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_sparse_matrix_multiplication/TestSparseMatrixSum/test_sparse_matrix_sum_small
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sparse_matrix_multiplication/TestSparseMatrixSum/test_sparse_matrix_sum_full
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __sparse_matrix_sum(self, file_name = "auxiliar_files_for_python_unittest/sparse_matrix_files/A.mm"): # Read the matrices`

# Calls

- [ReadMatrixMarketMatrix](../../../../functions/sources/matrix_market_interface/ReadMatrixMarketMatrix.md)
- [MatrixAdd](../../../../functions/utilities/sparse_matrix_multiplication_utility/MatrixAdd.md)

# Called by

- [test_sparse_matrix_sum_small](../../../../functions/tests/test_sparse_matrix_multiplication/TestSparseMatrixSum/test_sparse_matrix_sum_small.md)
- [test_sparse_matrix_sum_full](../../../../functions/tests/test_sparse_matrix_multiplication/TestSparseMatrixSum/test_sparse_matrix_sum_full.md)