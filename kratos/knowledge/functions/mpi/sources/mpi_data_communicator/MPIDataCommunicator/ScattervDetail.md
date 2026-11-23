---
type: C++ Method
title: ScattervDetail
resource: mpi/sources/mpi_data_communicator.cpp#L1026-L1061
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ValidateScattervInput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/CheckMPIErrorCode
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/includes/mpi_message/MPIMessage/SubDataTypeSize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/PrepareScattervBuffers
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MPIDataCommunicator::ScattervDetail( const TDataType& rSendValues, const std::vector<int>& rSendCounts, const std::vector<int>& rSendOffsets, TDataType& rRecvValues, const int SourceRank) const`

# Calls

- [ValidateScattervInput](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ValidateScattervInput.md)
- [CheckMPIErrorCode](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/CheckMPIErrorCode.md)
- [SubDataTypeSize](../../../../../functions/mpi/includes/mpi_message/MPIMessage/SubDataTypeSize.md)
- [PrepareScattervBuffers](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/PrepareScattervBuffers.md)