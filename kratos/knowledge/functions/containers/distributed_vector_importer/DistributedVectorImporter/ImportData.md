---
type: C++ Method
title: ImportData
resource: containers/distributed_vector_importer.h#L140-L173
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/containers/distributed_csr_matrix/DistributedCsrMatrix/ApplyHomogeneousDirichlet
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/distributed_csr_matrix/DistributedCsrMatrix/SpMV
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/cpp_tests/containers/test_distributed_system_vector/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_distributed_sparse_matrices/TestDistributedSparseMatrices/test_matrix_construction
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`DenseVector<TDataType> ImportData(const DistributedSystemVector<TDataType, TIndexType>& data_vector) const`

# Called by

- [ApplyHomogeneousDirichlet](../../../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/ApplyHomogeneousDirichlet.md)
- [SpMV](../../../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/SpMV.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/mpi/tests/cpp_tests/containers/test_distributed_system_vector/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_matrix_construction](../../../../functions/mpi/tests/test_distributed_sparse_matrices/TestDistributedSparseMatrices/test_matrix_construction.md)