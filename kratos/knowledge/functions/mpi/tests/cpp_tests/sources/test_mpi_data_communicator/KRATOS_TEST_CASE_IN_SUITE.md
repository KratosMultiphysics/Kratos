---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: mpi/tests/cpp_tests/sources/test_mpi_data_communicator.cpp#L28-L48
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/GetMPICommunicator
    resolved_by: tree-sitter
    confidence: exact
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

`KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorRankAndSize, KratosMPICoreFastSuite)`

# Calls

- [GetMPICommunicator](../../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/GetMPICommunicator.md)
- [Send](../../../../../../functions/includes/data_communicator/Send.md)
- [Recv](../../../../../../functions/includes/data_communicator/Recv.md)
- [Broadcast](../../../../../../functions/includes/data_communicator/Broadcast.md)