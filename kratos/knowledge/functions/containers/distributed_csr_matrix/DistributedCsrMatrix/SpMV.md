---
type: C++ Method
title: SpMV
resource: containers/distributed_csr_matrix.h#L493-L500
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/distributed_vector_importer/DistributedVectorImporter/ImportData
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SpMV(const DistributedSystemVector<TDataType,TIndexType>& x, DistributedSystemVector<TDataType,TIndexType>& y) const`

# Calls

- [ImportData](../../../../functions/containers/distributed_vector_importer/DistributedVectorImporter/ImportData.md)