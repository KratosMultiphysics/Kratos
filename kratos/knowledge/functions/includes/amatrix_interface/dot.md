---
type: C++ Function
title: dot
resource: includes/amatrix_interface.h#L250-L258
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/includes/amatrix_interface/Matrix/squared_norm
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/Matrix/norm
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_distributed_sparse_matrices/TestDistributedSparseMatrices/test_matrix_construction
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sparse_matrix_multiplication/TestSparseMatrixMultiplication/__sparse_matrix_multiplication
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`data_type dot( AMatrix::MatrixExpression<TExpressionType, AMatrix::row_major_access> const& Other) const`

# Called by

- [squared_norm](../../../functions/includes/amatrix_interface/Matrix/squared_norm.md)
- [norm](../../../functions/includes/amatrix_interface/Matrix/norm.md)
- [test_matrix_construction](../../../functions/mpi/tests/test_distributed_sparse_matrices/TestDistributedSparseMatrices/test_matrix_construction.md)
- [__sparse_matrix_multiplication](../../../functions/tests/test_sparse_matrix_multiplication/TestSparseMatrixMultiplication/__sparse_matrix_multiplication.md)