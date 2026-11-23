---
type: C++ Method
title: ScatterDetail
resource: mpi/sources/mpi_data_communicator.cpp#L970-L998
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/IsValidRank
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/IsEqualOnAllRanks
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/CheckMPIErrorCode
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/data_communicator/Broadcast
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MPIDataCommunicator::ScatterDetail( const TSendDataType& rSendValues, TRecvDataType& rRecvValues, const int SourceRank) const`

# Calls

- [IsValidRank](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/IsValidRank.md)
- [IsEqualOnAllRanks](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/IsEqualOnAllRanks.md)
- [CheckMPIErrorCode](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/CheckMPIErrorCode.md)
- [Broadcast](../../../../../functions/includes/data_communicator/Broadcast.md)