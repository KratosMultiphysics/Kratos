---
type: C++ Method
title: ApplyHomogeneousDirichlet
resource: containers/distributed_csr_matrix.h#L441-L478
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/distributed_csr_matrix/DistributedCsrMatrix/local_size1
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/distributed_vector_importer/DistributedVectorImporter/ImportData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/csr_matrix/CsrMatrix/nnz
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ApplyHomogeneousDirichlet(const DistributedSystemVector<TDataType,TIndexType>& rFreeDofsVector, const TDataType DiagonalValue, DistributedSystemVector<TDataType,TIndexType>& rRHS)`

# Calls

- [local_size1](../../../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/local_size1.md)
- [ImportData](../../../../functions/containers/distributed_vector_importer/DistributedVectorImporter/ImportData.md)
- [nnz](../../../../functions/containers/csr_matrix/CsrMatrix/nnz.md)