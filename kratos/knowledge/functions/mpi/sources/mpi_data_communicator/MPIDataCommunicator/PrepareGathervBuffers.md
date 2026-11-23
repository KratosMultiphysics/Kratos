---
type: C++ Method
title: PrepareGathervBuffers
resource: mpi/sources/mpi_data_communicator.cpp#L1521-L1555
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/GathervDetail
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MPIDataCommunicator::PrepareGathervBuffers( const std::vector<TDataType>& rGathervInput, std::vector<TDataType>& rGathervMessage, std::vector<int>& rMessageLengths, std::vector<int>& rMessageDistances, const int DestinationRank) const`

# Called by

- [GathervDetail](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/GathervDetail.md)