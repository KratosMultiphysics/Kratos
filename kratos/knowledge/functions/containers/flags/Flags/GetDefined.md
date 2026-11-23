---
type: C++ Method
title: GetDefined
resource: containers/flags.cpp#L84-L87
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/AndReduce
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/OrReduce
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/AndReduceAll
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/OrReduceAll
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`Flags::BlockType Flags::GetDefined() const`

# Called by

- [AndReduce](../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/AndReduce.md)
- [OrReduce](../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/OrReduce.md)
- [AndReduceAll](../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/AndReduceAll.md)
- [OrReduceAll](../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/OrReduceAll.md)