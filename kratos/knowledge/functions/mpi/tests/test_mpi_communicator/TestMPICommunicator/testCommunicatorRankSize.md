---
type: Python Method
title: testCommunicatorRankSize
resource: mpi/tests/test_mpi_communicator.py#L136-L149
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/mpi/tests/test_mpi_communicator/TestMPICommunicator/_read_model_part_mpi
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/MyPID
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/TotalProcesses
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def testCommunicatorRankSize(self):`

# Calls

- [_read_model_part_mpi](../../../../../functions/mpi/tests/test_mpi_communicator/TestMPICommunicator/_read_model_part_mpi.md)
- [MyPID](../../../../../functions/sources/communicator/Communicator/MyPID.md)
- [TotalProcesses](../../../../../functions/sources/communicator/Communicator/TotalProcesses.md)