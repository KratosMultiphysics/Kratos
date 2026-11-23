---
type: C++ Method
title: SubDataTypeSize
resource: mpi/includes/mpi_message.h#L193-L200
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ScattervDetail
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/GathervDetail
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/AllGathervDetail
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`inline int SubDataTypeSize(const MessageDataType& rValue)`

# Called by

- [ScattervDetail](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/ScattervDetail.md)
- [GathervDetail](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/GathervDetail.md)
- [AllGathervDetail](../../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/AllGathervDetail.md)