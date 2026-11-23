---
type: C++ Method
title: PrepareAllGathervBuffers
resource: mpi/sources/mpi_data_communicator.cpp#L1558-L1582
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/AllGathervDetail
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MPIDataCommunicator::PrepareAllGathervBuffers( const std::vector<TDataType>& rGathervInput, std::vector<TDataType>& rGathervMessage, std::vector<int>& rMessageLengths, std::vector<int>& rMessageDistances) const`

# Called by

- [AllGathervDetail](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/AllGathervDetail.md)