---
type: Python Method
title: __sparse_matrix_multiplication
resource: tests/test_sparse_matrix_multiplication.py#L76-L94
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/matrix_market_interface/ReadMatrixMarketMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/dot
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/sparse_matrix_multiplication_utility/MatrixMultiplicationSaad
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/sparse_matrix_multiplication_utility/MatrixMultiplicationRMerge
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_sparse_matrix_multiplication/TestSparseMatrixMultiplication/test_sparse_matrix_multiplication_saad_small
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sparse_matrix_multiplication/TestSparseMatrixMultiplication/test_sparse_matrix_multiplication_rmerge_small
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sparse_matrix_multiplication/TestSparseMatrixMultiplication/test_sparse_matrix_multiplication_saad_full
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sparse_matrix_multiplication/TestSparseMatrixMultiplication/test_sparse_matrix_multiplication_rmerge_full
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __sparse_matrix_multiplication(self, problem = "saad", file_name = "auxiliar_files_for_python_unittest/sparse_matrix_files/A.mm"): # Read the matrices`

# Calls

- [ReadMatrixMarketMatrix](../../../../functions/sources/matrix_market_interface/ReadMatrixMarketMatrix.md)
- [dot](../../../../functions/includes/amatrix_interface/dot.md)
- [MatrixMultiplicationSaad](../../../../functions/utilities/sparse_matrix_multiplication_utility/MatrixMultiplicationSaad.md)
- [MatrixMultiplicationRMerge](../../../../functions/utilities/sparse_matrix_multiplication_utility/MatrixMultiplicationRMerge.md)

# Called by

- [test_sparse_matrix_multiplication_saad_small](../../../../functions/tests/test_sparse_matrix_multiplication/TestSparseMatrixMultiplication/test_sparse_matrix_multiplication_saad_small.md)
- [test_sparse_matrix_multiplication_rmerge_small](../../../../functions/tests/test_sparse_matrix_multiplication/TestSparseMatrixMultiplication/test_sparse_matrix_multiplication_rmerge_small.md)
- [test_sparse_matrix_multiplication_saad_full](../../../../functions/tests/test_sparse_matrix_multiplication/TestSparseMatrixMultiplication/test_sparse_matrix_multiplication_saad_full.md)
- [test_sparse_matrix_multiplication_rmerge_full](../../../../functions/tests/test_sparse_matrix_multiplication/TestSparseMatrixMultiplication/test_sparse_matrix_multiplication_rmerge_full.md)