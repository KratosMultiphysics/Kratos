---
type: C++ Method
title: RemoteLocalId
resource: containers/distributed_numbering.h#L189-L192
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/containers/distributed_sparse_graph/DistributedSparseGraph/AddEntry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/distributed_sparse_graph/AddEntries
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/distributed_vector_exporter/DistributedVectorExporter
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/distributed_vector_importer/DistributedVectorImporter
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`inline IndexType RemoteLocalId(const IndexType rGlobalId, const IndexType rOwnerRank) const`

# Called by

- [AddEntry](../../../../functions/containers/distributed_sparse_graph/DistributedSparseGraph/AddEntry.md)
- [AddEntries](../../../../functions/containers/distributed_sparse_graph/AddEntries.md)
- [DistributedVectorExporter](../../../../functions/containers/distributed_vector_exporter/DistributedVectorExporter.md)
- [DistributedVectorImporter](../../../../functions/containers/distributed_vector_importer/DistributedVectorImporter.md)