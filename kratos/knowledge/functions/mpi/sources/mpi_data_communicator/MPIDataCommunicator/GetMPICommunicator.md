---
type: C++ Method
title: GetMPICommunicator
resource: mpi/sources/mpi_data_communicator.cpp#L431-L442
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/mpi/includes/mpi_communicator/AsyncSendAndReceiveObjects
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/cpp_tests/sources/test_mpi_data_communicator/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/amgcl_distributed_csr_conversion_utilities/ConvertToAmgcl
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/amgcl_distributed_csr_conversion_utilities/ConvertToCsrMatrix
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`MPI_Comm MPIDataCommunicator::GetMPICommunicator(const DataCommunicator& rDataCommunicator)`

# Called by

- [AsyncSendAndReceiveObjects](../../../../../functions/mpi/includes/mpi_communicator/AsyncSendAndReceiveObjects.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/mpi/tests/cpp_tests/sources/test_mpi_data_communicator/KRATOS_TEST_CASE_IN_SUITE.md)
- [ConvertToAmgcl](../../../../../functions/mpi/utilities/amgcl_distributed_csr_conversion_utilities/ConvertToAmgcl.md)
- [ConvertToCsrMatrix](../../../../../functions/mpi/utilities/amgcl_distributed_csr_conversion_utilities/ConvertToCsrMatrix.md)