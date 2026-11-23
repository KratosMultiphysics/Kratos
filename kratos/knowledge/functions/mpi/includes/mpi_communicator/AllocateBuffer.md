---
type: C++ Function
title: AllocateBuffer
resource: mpi/includes/mpi_communicator.h#L1917-L1925
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/mpi/includes/mpi_communicator/TransferDistributedValues
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AllocateBuffer(std::vector<TSendType>& rBuffer, const MeshType& rSourceMesh, TDatabaseAccess& rAccess)`

# Called by

- [TransferDistributedValues](../../../../functions/mpi/includes/mpi_communicator/TransferDistributedValues.md)