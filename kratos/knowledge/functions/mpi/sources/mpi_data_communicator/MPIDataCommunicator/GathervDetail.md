---
type: C++ Method
title: GathervDetail
resource: mpi/sources/mpi_data_communicator.cpp#L1128-L1167
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ValidateGathervInput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/CheckMPIErrorCode
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/includes/mpi_message/MPIMessage/SubDataTypeSize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/PrepareGathervBuffers
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/PrepareGathervReturn
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MPIDataCommunicator::GathervDetail( const TDataType& rSendValues, TDataType& rRecvValues, const std::vector<int>& rRecvCounts, const std::vector<int>& rRecvOffsets, const int RecvRank) const`

# Calls

- [ValidateGathervInput](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ValidateGathervInput.md)
- [CheckMPIErrorCode](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/CheckMPIErrorCode.md)
- [SubDataTypeSize](../../../../../functions/mpi/includes/mpi_message/MPIMessage/SubDataTypeSize.md)
- [PrepareGathervBuffers](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/PrepareGathervBuffers.md)
- [PrepareGathervReturn](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/PrepareGathervReturn.md)