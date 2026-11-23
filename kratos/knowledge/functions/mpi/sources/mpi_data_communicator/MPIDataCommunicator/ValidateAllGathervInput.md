---
type: C++ Method
title: ValidateAllGathervInput
resource: mpi/sources/mpi_data_communicator.cpp#L1433-L1472
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/CheckMPIErrorCode
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/AllGathervDetail
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MPIDataCommunicator::ValidateAllGathervInput( const TDataType& rSendValues, TDataType& rRecvValues, const std::vector<int>& rRecvCounts, const std::vector<int>& rRecvOffsets) const`

# Calls

- [CheckMPIErrorCode](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/CheckMPIErrorCode.md)

# Called by

- [AllGathervDetail](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/AllGathervDetail.md)