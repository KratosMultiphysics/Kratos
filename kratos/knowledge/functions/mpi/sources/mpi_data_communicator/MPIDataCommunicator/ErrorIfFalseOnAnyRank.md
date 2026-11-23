---
type: C++ Method
title: ErrorIfFalseOnAnyRank
resource: mpi/sources/mpi_data_communicator.cpp#L1311-L1319
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/CheckMPIErrorCode
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool MPIDataCommunicator::ErrorIfFalseOnAnyRank(bool Condition) const`

# Calls

- [CheckMPIErrorCode](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/CheckMPIErrorCode.md)