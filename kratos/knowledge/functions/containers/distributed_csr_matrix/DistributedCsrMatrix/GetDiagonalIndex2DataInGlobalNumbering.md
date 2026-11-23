---
type: C++ Method
title: GetDiagonalIndex2DataInGlobalNumbering
resource: containers/distributed_csr_matrix.h#L421-L429
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/distributed_numbering/DistributedNumbering/GlobalId
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/containers/distributed_csr_matrix/DistributedCsrMatrix/PrintData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/python/add_distributed_sparse_matrices_to_python/AddDistributedSparseMatricesToPython
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`DenseVector<TIndexType> GetDiagonalIndex2DataInGlobalNumbering() const`

# Calls

- [GlobalId](../../../../functions/containers/distributed_numbering/DistributedNumbering/GlobalId.md)

# Called by

- [PrintData](../../../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/PrintData.md)
- [AddDistributedSparseMatricesToPython](../../../../functions/mpi/python/add_distributed_sparse_matrices_to_python/AddDistributedSparseMatricesToPython.md)