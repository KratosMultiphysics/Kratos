---
type: C++ Method
title: GlobalId
resource: containers/distributed_numbering.h#L183-L187
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/containers/distributed_csr_matrix/DistributedCsrMatrix/GetDiagonalIndex2DataInGlobalNumbering
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/distributed_csr_matrix/DistributedCsrMatrix/ToMap
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/distributed_csr_matrix/DistributedCsrMatrix/ToSerialCSR
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/cpp_tests/containers/test_distributed_csr_matrix/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/cpp_tests/containers/test_distributed_system_vector/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_distributed_sparse_matrices/TestDistributedSparseMatrices/test_matrix_construction
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_utilities/distributed_sparse_containers_test_utilities/CheckGraph
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`inline IndexType GlobalId(const IndexType rLocalId) const`

# Called by

- [GetDiagonalIndex2DataInGlobalNumbering](../../../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/GetDiagonalIndex2DataInGlobalNumbering.md)
- [ToMap](../../../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/ToMap.md)
- [ToSerialCSR](../../../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/ToSerialCSR.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/mpi/tests/cpp_tests/containers/test_distributed_csr_matrix/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/mpi/tests/cpp_tests/containers/test_distributed_system_vector/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_matrix_construction](../../../../functions/mpi/tests/test_distributed_sparse_matrices/TestDistributedSparseMatrices/test_matrix_construction.md)
- [CheckGraph](../../../../functions/mpi/tests/test_utilities/distributed_sparse_containers_test_utilities/CheckGraph.md)