---
type: C++ Function
title: AsyncSendAndReceiveObjects
resource: mpi/includes/mpi_communicator.h#L1448-L1560
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/GetMPICommunicator
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/serializer/pGetBuffer
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool AsyncSendAndReceiveObjects(std::vector<TObjectType>& SendObjects, std::vector<TObjectType>& RecvObjects)`

# Calls

- [GetMPICommunicator](../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/GetMPICommunicator.md)
- [pGetBuffer](../../../../functions/includes/serializer/pGetBuffer.md)