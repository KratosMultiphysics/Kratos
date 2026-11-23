---
type: C++ Method
title: PrepareScattervBuffers
resource: mpi/sources/mpi_data_communicator.cpp#L1474-L1518
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ScattervDetail
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MPIDataCommunicator::PrepareScattervBuffers( const std::vector<std::vector<TDataType>>& rInputMessage, std::vector<TDataType>& rScattervMessage, std::vector<int>& rMessageLengths, std::vector<int>& rMessageDistances, std::vector<TDataType>& rResult, const int SourceRank) const`

# Called by

- [ScattervDetail](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ScattervDetail.md)