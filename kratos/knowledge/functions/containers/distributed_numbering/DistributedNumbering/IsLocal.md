---
type: C++ Method
title: IsLocal
resource: containers/distributed_numbering.h#L171-L175
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/containers/distributed_csr_matrix/DistributedCsrMatrix/DistributedCsrMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/distributed_csr_matrix/Assemble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/distributed_csr_matrix/DistributedCsrMatrix/AssembleEntry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/distributed_sparse_graph/DistributedSparseGraph/AddEntry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/distributed_sparse_graph/AddEntries
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/distributed_system_vector/DistributedSystemVector/AddEntry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/distributed_system_vector/Assemble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/cpp_tests/containers/test_distributed_csr_matrix/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/cpp_tests/containers/test_distributed_exporter/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`inline bool IsLocal(const IndexType I) const`

# Called by

- [DistributedCsrMatrix](../../../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/DistributedCsrMatrix.md)
- [Assemble](../../../../functions/containers/distributed_csr_matrix/Assemble.md)
- [AssembleEntry](../../../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/AssembleEntry.md)
- [AddEntry](../../../../functions/containers/distributed_sparse_graph/DistributedSparseGraph/AddEntry.md)
- [AddEntries](../../../../functions/containers/distributed_sparse_graph/AddEntries.md)
- [AddEntry](../../../../functions/containers/distributed_system_vector/DistributedSystemVector/AddEntry.md)
- [Assemble](../../../../functions/containers/distributed_system_vector/Assemble.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/mpi/tests/cpp_tests/containers/test_distributed_csr_matrix/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/mpi/tests/cpp_tests/containers/test_distributed_exporter/KRATOS_TEST_CASE_IN_SUITE.md)