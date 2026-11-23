---
type: C++ Function
title: ConvertToAmgcl
resource: mpi/utilities/amgcl_distributed_csr_conversion_utilities.h#L47-L77
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/distributed_csr_matrix/DistributedCsrMatrix/local_size1
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/GetMPICommunicator
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static Kratos::shared_ptr<amgcl::mpi::distributed_matrix<amgcl::backend::builtin<double>>> ConvertToAmgcl( const DistributedCsrMatrix<TDataType, TIndexType>& rA, DenseVector<TIndexType>& global_index2, bool MoveToBackend=true)`

# Calls

- [local_size1](../../../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/local_size1.md)
- [GetMPICommunicator](../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/GetMPICommunicator.md)