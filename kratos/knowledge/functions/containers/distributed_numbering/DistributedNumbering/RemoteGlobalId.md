---
type: C++ Method
title: RemoteGlobalId
resource: containers/distributed_numbering.h#L194-L197
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/containers/distributed_csr_matrix/DistributedCsrMatrix/PrepareNonLocalCommunications
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/distributed_system_vector/DistributedSystemVector/DistributedSystemVector
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`inline IndexType RemoteGlobalId(const IndexType rRemoteLocalId, const IndexType rOwnerRank) const`

# Called by

- [PrepareNonLocalCommunications](../../../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/PrepareNonLocalCommunications.md)
- [DistributedSystemVector](../../../../functions/containers/distributed_system_vector/DistributedSystemVector/DistributedSystemVector.md)