---
type: C++ Method
title: CheckMPIErrorCode
resource: mpi/sources/mpi_data_communicator.cpp#L524-L527
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/Barrier
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/Rank
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/Size
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/SynchronizeShapeDetail
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ReduceDetail
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/AllReduceDetail
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ScanDetail
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/SendRecvDetail
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/SendDetail
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/RecvDetail
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/BroadcastDetail
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ScatterDetail
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ScattervDetail
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/GatherDetail
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/GathervDetail
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/AllGatherDetail
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/AllGathervDetail
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/BroadcastErrorIfTrue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/BroadcastErrorIfFalse
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ErrorIfTrueOnAnyRank
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ErrorIfFalseOnAnyRank
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/IsEqualOnAllRanks
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ValidateScattervInput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ValidateGathervInput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ValidateAllGathervInput
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MPIDataCommunicator::CheckMPIErrorCode(const int ierr, const std::string& MPICallName) const`

# Called by

- [Barrier](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/Barrier.md)
- [Rank](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/Rank.md)
- [Size](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/Size.md)
- [SynchronizeShapeDetail](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/SynchronizeShapeDetail.md)
- [ReduceDetail](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ReduceDetail.md)
- [AllReduceDetail](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/AllReduceDetail.md)
- [ScanDetail](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ScanDetail.md)
- [SendRecvDetail](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/SendRecvDetail.md)
- [SendDetail](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/SendDetail.md)
- [RecvDetail](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/RecvDetail.md)
- [BroadcastDetail](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/BroadcastDetail.md)
- [ScatterDetail](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ScatterDetail.md)
- [ScattervDetail](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ScattervDetail.md)
- [GatherDetail](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/GatherDetail.md)
- [GathervDetail](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/GathervDetail.md)
- [AllGatherDetail](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/AllGatherDetail.md)
- [AllGathervDetail](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/AllGathervDetail.md)
- [BroadcastErrorIfTrue](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/BroadcastErrorIfTrue.md)
- [BroadcastErrorIfFalse](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/BroadcastErrorIfFalse.md)
- [ErrorIfTrueOnAnyRank](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ErrorIfTrueOnAnyRank.md)
- [ErrorIfFalseOnAnyRank](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ErrorIfFalseOnAnyRank.md)
- [IsEqualOnAllRanks](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/IsEqualOnAllRanks.md)
- [ValidateScattervInput](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ValidateScattervInput.md)
- [ValidateGathervInput](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ValidateGathervInput.md)
- [ValidateAllGathervInput](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ValidateAllGathervInput.md)