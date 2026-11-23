---
type: C++ Method
title: PrepareNonLocalCommunications
resource: containers/distributed_csr_matrix.h#L922-L966
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/communication_coloring_utilities/MPIColoringUtilities/ComputeCommunicationScheduling
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/distributed_numbering/DistributedNumbering/RemoteGlobalId
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/containers/distributed_csr_matrix/DistributedCsrMatrix/DistributedCsrMatrix
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void PrepareNonLocalCommunications(const DistributedSparseGraph<TIndexType>& rSparseGraph)`

# Calls

- [ComputeCommunicationScheduling](../../../../functions/utilities/communication_coloring_utilities/MPIColoringUtilities/ComputeCommunicationScheduling.md)
- [RemoteGlobalId](../../../../functions/containers/distributed_numbering/DistributedNumbering/RemoteGlobalId.md)

# Called by

- [DistributedCsrMatrix](../../../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/DistributedCsrMatrix.md)