---
type: C++ Function
title: ConvertToCsrMatrix
resource: mpi/utilities/amgcl_distributed_csr_conversion_utilities.h#L80-L141
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/GetMPICommunicator
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/csr_matrix/CsrMatrix/SetColSize
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static typename DistributedCsrMatrix<TDataType, TIndexType>::UniquePointer ConvertToCsrMatrix( amgcl::mpi::distributed_matrix<amgcl::backend::builtin<double>>& rA, //cannot be made const since i need to modify some data in-place, const DataCommunicator& kratos_comm )`

# Calls

- [GetMPICommunicator](../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/GetMPICommunicator.md)
- [SetColSize](../../../../functions/containers/csr_matrix/CsrMatrix/SetColSize.md)