---
type: C++ Method
title: AllGatherDetail
resource: mpi/sources/mpi_data_communicator.cpp#L1184-L1209
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/IsEqualOnAllRanks
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/CheckMPIErrorCode
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MPIDataCommunicator::AllGatherDetail( const TDataType& rSendValues, TDataType& rRecvValues) const`

# Calls

- [IsEqualOnAllRanks](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/IsEqualOnAllRanks.md)
- [CheckMPIErrorCode](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/CheckMPIErrorCode.md)