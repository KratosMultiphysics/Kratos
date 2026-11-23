---
type: C++ Method
title: PrepareGathervReturn
resource: mpi/sources/mpi_data_communicator.cpp#L1584-L1604
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

`void MPIDataCommunicator::PrepareGathervReturn( const std::vector<TDataType>& rGathervMessage, const std::vector<int>& rMessageLengths, const std::vector<int>& rMessageDistances, std::vector<std::vector<TDataType>>& rOutputMessage, const int DestinationRank) const`

# Called by

- [GathervDetail](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/GathervDetail.md)