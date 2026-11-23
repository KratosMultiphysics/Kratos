---
type: C++ Method
title: GetStringRepresentation
resource: includes/stream_serializer.h#L53-L55
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/serializer/pGetBuffer
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/includes/data_communicator/BroadcastImpl
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/data_communicator/SendRecvImpl
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/data_communicator/SendImpl
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/includes/mpi_communicator/FillBuffer
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_serializer_to_python/AddSerializerToPython
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::string GetStringRepresentation()`

# Calls

- [pGetBuffer](../../../../functions/includes/serializer/pGetBuffer.md)

# Called by

- [BroadcastImpl](../../../../functions/includes/data_communicator/BroadcastImpl.md)
- [SendRecvImpl](../../../../functions/includes/data_communicator/SendRecvImpl.md)
- [SendImpl](../../../../functions/includes/data_communicator/SendImpl.md)
- [FillBuffer](../../../../functions/mpi/includes/mpi_communicator/FillBuffer.md)
- [AddSerializerToPython](../../../../functions/python/add_serializer_to_python/AddSerializerToPython.md)