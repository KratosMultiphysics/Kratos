---
type: C++ Method
title: ReduceDetail
resource: mpi/sources/mpi_data_communicator.cpp#L632-L661
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
---

# Signature

`void MPIDataCommunicator::ReduceDetail( const TDataType& rLocalValues, TDataType& rReducedValues, MPI_Op Operation, const int Root) const`

# Calls

- [IsValidRank](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/IsValidRank.md)
- [IsEqualOnAllRanks](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/IsEqualOnAllRanks.md)
- [CheckMPIErrorCode](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/CheckMPIErrorCode.md)