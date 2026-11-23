---
type: C++ Method
title: DistributedSystemVector
resource: containers/distributed_system_vector.h#L76-L96
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/distributed_numbering/DistributedNumbering/RemoteGlobalId
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`DistributedSystemVector(const DistributedSparseGraph<IndexType>& rGraph) : mpComm(&rGraph.GetComm())`

# Calls

- [RemoteGlobalId](../../../../functions/containers/distributed_numbering/DistributedNumbering/RemoteGlobalId.md)