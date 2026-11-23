---
type: C++ Method
title: AddEntry
resource: containers/distributed_sparse_graph.h#L168-L180
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/distributed_numbering/DistributedNumbering/IsLocal
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/distributed_numbering/DistributedNumbering/LocalId
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/distributed_numbering/DistributedNumbering/OwnerRank
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/distributed_numbering/DistributedNumbering/RemoteLocalId
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AddEntry(const IndexType RowIndex, const IndexType ColIndex)`

# Calls

- [IsLocal](../../../../functions/containers/distributed_numbering/DistributedNumbering/IsLocal.md)
- [LocalId](../../../../functions/containers/distributed_numbering/DistributedNumbering/LocalId.md)
- [OwnerRank](../../../../functions/containers/distributed_numbering/DistributedNumbering/OwnerRank.md)
- [RemoteLocalId](../../../../functions/containers/distributed_numbering/DistributedNumbering/RemoteLocalId.md)