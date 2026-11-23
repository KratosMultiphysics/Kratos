---
type: C++ Method
title: BroadcastErrorIfTrue
resource: mpi/sources/mpi_data_communicator.cpp#L1277-L1285
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

`bool MPIDataCommunicator::BroadcastErrorIfTrue(bool Condition, const int SourceRank) const`

# Calls

- [CheckMPIErrorCode](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/CheckMPIErrorCode.md)