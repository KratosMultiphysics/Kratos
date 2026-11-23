---
type: C++ Method
title: ToSerialCSR
resource: containers/distributed_csr_matrix.h#L700-L766
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/csr_matrix/CsrMatrix/nnz
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/distributed_numbering/DistributedNumbering/GlobalId
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/tests/test_distributed_sparse_matrices/TestDistributedSparseMatrices/test_matrix_construction
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`typename CsrMatrix<TDataType,TIndexType>::Pointer ToSerialCSR(MpiIndexType target_rank=0) const`

# Calls

- [nnz](../../../../functions/containers/csr_matrix/CsrMatrix/nnz.md)
- [GlobalId](../../../../functions/containers/distributed_numbering/DistributedNumbering/GlobalId.md)

# Called by

- [test_matrix_construction](../../../../functions/mpi/tests/test_distributed_sparse_matrices/TestDistributedSparseMatrices/test_matrix_construction.md)