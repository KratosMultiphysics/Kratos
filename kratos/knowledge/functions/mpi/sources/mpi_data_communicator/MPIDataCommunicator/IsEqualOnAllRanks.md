---
type: C++ Method
title: IsEqualOnAllRanks
resource: mpi/sources/mpi_data_communicator.cpp#L1321-L1330
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/CheckMPIErrorCode
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ReduceDetail
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/AllReduceDetail
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ScanDetail
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/BroadcastDetail
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ScatterDetail
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/GatherDetail
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/AllGatherDetail
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool MPIDataCommunicator::IsEqualOnAllRanks(const int LocalValue) const`

# Calls

- [CheckMPIErrorCode](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/CheckMPIErrorCode.md)

# Called by

- [ReduceDetail](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ReduceDetail.md)
- [AllReduceDetail](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/AllReduceDetail.md)
- [ScanDetail](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ScanDetail.md)
- [BroadcastDetail](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/BroadcastDetail.md)
- [ScatterDetail](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ScatterDetail.md)
- [GatherDetail](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/GatherDetail.md)
- [AllGatherDetail](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/AllGatherDetail.md)