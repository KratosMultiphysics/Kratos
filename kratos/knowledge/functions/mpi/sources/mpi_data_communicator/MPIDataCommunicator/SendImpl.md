---
type: C++ Method
title: SendImpl
resource: mpi/sources/mpi_data_communicator.cpp#L580-L583
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/SendDetail
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MPIDataCommunicator::SendImpl(const std::string& rSendValues, const int SendDestination, const int SendTag) const`

# Calls

- [SendDetail](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/SendDetail.md)