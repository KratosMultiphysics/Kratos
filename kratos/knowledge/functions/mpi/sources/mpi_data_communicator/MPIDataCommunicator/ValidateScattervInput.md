---
type: C++ Method
title: ValidateScattervInput
resource: mpi/sources/mpi_data_communicator.cpp#L1337-L1382
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/IsValidRank
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/CheckMPIErrorCode
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ScattervDetail
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MPIDataCommunicator::ValidateScattervInput( const TDataType& rSendValues, const std::vector<int>& rSendCounts, const std::vector<int>& rSendOffsets, TDataType& rRecvValues, const int SourceRank) const`

# Calls

- [IsValidRank](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/IsValidRank.md)
- [CheckMPIErrorCode](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/CheckMPIErrorCode.md)

# Called by

- [ScattervDetail](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ScattervDetail.md)