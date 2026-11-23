---
type: C++ Function
title: ReadMatrixMarketMatrix
resource: sources/matrix_market_interface.cpp#L82-L302
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/mmio/mm_read_banner
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/mmio/mm_is_valid
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/mmio/mm_typecode_to_str
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/mmio/mm_read_mtx_crd_size
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_condition_number/TestConditionNumber/test_condition_number
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_eigen_solvers/TestEigenSolvers/_auxiliary_test_function
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_linear_solvers/TestLinearSolvers/_auxiliary_test_function
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_matrix_market_interface/TestMatrixMarketInterface/test_mm_matrix_io
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_matrix_market_interface/TestMatrixMarketInterface/test_mm_matrix_io_cplx
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sparse_matrix_multiplication/TestSparseMatrixSum/__sparse_matrix_sum
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sparse_matrix_multiplication/TestSparseMatrixTranspose/__sparse_matrix_transpose
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sparse_matrix_multiplication/TestSparseMatrixMultiplication/__sparse_matrix_multiplication
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool ReadMatrixMarketMatrix(const char *FileName, CompressedMatrixType &M)`

# Calls

- [mm_read_banner](../../../functions/sources/mmio/mm_read_banner.md)
- [mm_is_valid](../../../functions/sources/mmio/mm_is_valid.md)
- [mm_typecode_to_str](../../../functions/sources/mmio/mm_typecode_to_str.md)
- [mm_read_mtx_crd_size](../../../functions/sources/mmio/mm_read_mtx_crd_size.md)

# Called by

- [test_condition_number](../../../functions/tests/test_condition_number/TestConditionNumber/test_condition_number.md)
- [_auxiliary_test_function](../../../functions/tests/test_eigen_solvers/TestEigenSolvers/_auxiliary_test_function.md)
- [_auxiliary_test_function](../../../functions/tests/test_linear_solvers/TestLinearSolvers/_auxiliary_test_function.md)
- [test_mm_matrix_io](../../../functions/tests/test_matrix_market_interface/TestMatrixMarketInterface/test_mm_matrix_io.md)
- [test_mm_matrix_io_cplx](../../../functions/tests/test_matrix_market_interface/TestMatrixMarketInterface/test_mm_matrix_io_cplx.md)
- [__sparse_matrix_sum](../../../functions/tests/test_sparse_matrix_multiplication/TestSparseMatrixSum/__sparse_matrix_sum.md)
- [__sparse_matrix_transpose](../../../functions/tests/test_sparse_matrix_multiplication/TestSparseMatrixTranspose/__sparse_matrix_transpose.md)
- [__sparse_matrix_multiplication](../../../functions/tests/test_sparse_matrix_multiplication/TestSparseMatrixMultiplication/__sparse_matrix_multiplication.md)