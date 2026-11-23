---
type: C++ Method
title: BroadcastErrorIfFalse
resource: mpi/sources/mpi_data_communicator.cpp#L1287-L1295
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

`bool MPIDataCommunicator::BroadcastErrorIfFalse(bool Condition, const int SourceRank) const`

# Calls

- [CheckMPIErrorCode](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/CheckMPIErrorCode.md)