---
type: C++ Method
title: PrepareAllGathervReturn
resource: mpi/sources/mpi_data_communicator.cpp#L1607-L1621
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

`void MPIDataCommunicator::PrepareAllGathervReturn( const std::vector<TDataType>& rGathervMessage, const std::vector<int>& rMessageLengths, const std::vector<int>& rMessageDistances, std::vector<std::vector<TDataType>>& rOutputMessage) const`

# Called by

- [AllGathervDetail](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/AllGathervDetail.md)