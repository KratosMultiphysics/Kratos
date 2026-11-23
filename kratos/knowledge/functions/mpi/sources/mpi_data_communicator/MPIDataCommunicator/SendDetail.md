---
type: C++ Method
title: SendDetail
resource: mpi/sources/mpi_data_communicator.cpp#L870-L895
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/CheckMPIErrorCode
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/SendImpl
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MPIDataCommunicator::SendDetail( const TDataType& rSendValues, const int SendDestination, const int SendTag) const`

# Calls

- [CheckMPIErrorCode](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/CheckMPIErrorCode.md)

# Called by

- [SendImpl](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/SendImpl.md)