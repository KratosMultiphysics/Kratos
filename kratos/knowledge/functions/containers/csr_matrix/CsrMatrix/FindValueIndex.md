---
type: C++ Method
title: FindValueIndex
resource: containers/csr_matrix.h#L400-L405
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/containers/csr_matrix/CsrMatrix/Has
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_sparse_matrices_to_python/GetEquationIdCsrIndices
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/containers/test_csr_matrix/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`IndexType FindValueIndex(IndexType I, IndexType J) const`

# Called by

- [Has](../../../../functions/containers/csr_matrix/CsrMatrix/Has.md)
- [GetEquationIdCsrIndices](../../../../functions/python/add_sparse_matrices_to_python/GetEquationIdCsrIndices.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/containers/test_csr_matrix/KRATOS_TEST_CASE_IN_SUITE.md)