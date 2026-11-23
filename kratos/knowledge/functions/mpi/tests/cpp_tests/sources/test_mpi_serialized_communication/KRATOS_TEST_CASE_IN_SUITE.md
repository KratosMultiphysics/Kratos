---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: mpi/tests/cpp_tests/sources/test_mpi_serialized_communication.cpp#L26-L51
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/data_communicator/Send
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/data_communicator/Recv
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/data_communicator/Broadcast
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorSerializedSendRecv, KratosMPICoreFastSuite)`

# Calls

- [Send](../../../../../../functions/includes/data_communicator/Send.md)
- [Recv](../../../../../../functions/includes/data_communicator/Recv.md)
- [Broadcast](../../../../../../functions/includes/data_communicator/Broadcast.md)