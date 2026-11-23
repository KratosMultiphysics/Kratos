---
type: C++ Method
title: AllGathervDetail
resource: mpi/sources/mpi_data_communicator.cpp#L1225-L1261
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ValidateAllGathervInput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/CheckMPIErrorCode
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/includes/mpi_message/MPIMessage/SubDataTypeSize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/PrepareAllGathervBuffers
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/PrepareAllGathervReturn
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MPIDataCommunicator::AllGathervDetail( const TDataType& rSendValues, TDataType& rRecvValues, const std::vector<int>& rRecvCounts, const std::vector<int>& rRecvOffsets) const`

# Calls

- [ValidateAllGathervInput](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ValidateAllGathervInput.md)
- [CheckMPIErrorCode](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/CheckMPIErrorCode.md)
- [SubDataTypeSize](../../../../../functions/mpi/includes/mpi_message/MPIMessage/SubDataTypeSize.md)
- [PrepareAllGathervBuffers](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/PrepareAllGathervBuffers.md)
- [PrepareAllGathervReturn](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/PrepareAllGathervReturn.md)