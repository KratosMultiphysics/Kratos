---
type: C++ Method
title: RecvImpl
resource: mpi/sources/mpi_data_communicator.cpp#L585-L588
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/RecvDetail
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MPIDataCommunicator::RecvImpl(std::string& rRecvValues, const int RecvSource, const int RecvTag) const`

# Calls

- [RecvDetail](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/RecvDetail.md)