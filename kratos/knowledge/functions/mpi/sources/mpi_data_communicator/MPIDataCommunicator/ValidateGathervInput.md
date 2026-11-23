---
type: C++ Method
title: ValidateGathervInput
resource: mpi/sources/mpi_data_communicator.cpp#L1384-L1430
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
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/GathervDetail
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MPIDataCommunicator::ValidateGathervInput( const TDataType& rSendValues, TDataType& rRecvValues, const std::vector<int>& rRecvCounts, const std::vector<int>& rRecvOffsets, const int RecvRank) const`

# Calls

- [IsValidRank](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/IsValidRank.md)
- [CheckMPIErrorCode](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/CheckMPIErrorCode.md)

# Called by

- [GathervDetail](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/GathervDetail.md)