---
type: C++ Method
title: SendRecvDetail
resource: mpi/sources/mpi_data_communicator.cpp#L816-L829
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/CheckMPIErrorCode
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MPIDataCommunicator::SendRecvDetail( const TDataType& rSendMessage, const int SendDestination, const int SendTag, TDataType& rRecvMessage, const int RecvSource, const int RecvTag) const`

# Calls

- [CheckMPIErrorCode](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/CheckMPIErrorCode.md)