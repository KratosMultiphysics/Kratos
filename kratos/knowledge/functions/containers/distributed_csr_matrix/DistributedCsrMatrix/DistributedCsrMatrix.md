---
type: C++ Method
title: DistributedCsrMatrix
resource: containers/distributed_csr_matrix.h#L74-L75
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/distributed_sparse_graph/DistributedSparseGraph/ComputeMaxGlobalColumnIndex
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/distributed_numbering/DistributedNumbering/IsLocal
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/csr_matrix/CsrMatrix/SetColSize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/distributed_numbering/DistributedNumbering/LocalId
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/csr_matrix/CsrMatrix/CheckColSize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/distributed_csr_matrix/DistributedCsrMatrix/GetOffDiagonalBlockLocalId
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/distributed_csr_matrix/DistributedCsrMatrix/PrepareNonLocalCommunications
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/distributed_csr_matrix/DistributedCsrMatrix/ReconstructDirectAccessVectors
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`DistributedCsrMatrix()`

# Calls

- [ComputeMaxGlobalColumnIndex](../../../../functions/containers/distributed_sparse_graph/DistributedSparseGraph/ComputeMaxGlobalColumnIndex.md)
- [IsLocal](../../../../functions/containers/distributed_numbering/DistributedNumbering/IsLocal.md)
- [SetColSize](../../../../functions/containers/csr_matrix/CsrMatrix/SetColSize.md)
- [LocalId](../../../../functions/containers/distributed_numbering/DistributedNumbering/LocalId.md)
- [CheckColSize](../../../../functions/containers/csr_matrix/CsrMatrix/CheckColSize.md)
- [GetOffDiagonalBlockLocalId](../../../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/GetOffDiagonalBlockLocalId.md)
- [PrepareNonLocalCommunications](../../../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/PrepareNonLocalCommunications.md)
- [ReconstructDirectAccessVectors](../../../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/ReconstructDirectAccessVectors.md)