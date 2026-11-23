---
type: C++ Method
title: IsValidRank
resource: mpi/sources/mpi_data_communicator.cpp#L1332-L1335
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ReduceDetail
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
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ValidateScattervInput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ValidateGathervInput
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool MPIDataCommunicator::IsValidRank(const int Rank) const`

# Called by

- [ReduceDetail](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ReduceDetail.md)
- [BroadcastDetail](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/BroadcastDetail.md)
- [ScatterDetail](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ScatterDetail.md)
- [GatherDetail](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/GatherDetail.md)
- [ValidateScattervInput](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ValidateScattervInput.md)
- [ValidateGathervInput](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ValidateGathervInput.md)