---
type: C++ Function
title: DistributedVectorImporter
resource: containers/distributed_vector_importer.h#L73-L126
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/distributed_numbering/DistributedNumbering/OwnerRank
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/distributed_numbering/DistributedNumbering/RemoteLocalId
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/communication_coloring_utilities/MPIColoringUtilities/ComputeCommunicationScheduling
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`DistributedVectorImporter( const DataCommunicator& rComm, const TGlobalIndicesVectorType& rGlobalIndices, const DistributedNumbering<IndexType>& rNumbering) : mrComm(rComm)`

# Calls

- [OwnerRank](../../../functions/containers/distributed_numbering/DistributedNumbering/OwnerRank.md)
- [RemoteLocalId](../../../functions/containers/distributed_numbering/DistributedNumbering/RemoteLocalId.md)
- [ComputeCommunicationScheduling](../../../functions/utilities/communication_coloring_utilities/MPIColoringUtilities/ComputeCommunicationScheduling.md)