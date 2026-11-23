---
type: C++ Method
title: CheckCSRMatrix
resource: tests/test_utilities/sparse_containers_test_utilities.cpp#L181-L211
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/csr_matrix/CsrMatrix/nnz
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/containers/test_csr_matrix/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool SparseContainersTestUtilities::CheckCSRMatrix( const CsrMatrix<>& rA, const MatrixMapType& rReferenceGraph)`

# Calls

- [nnz](../../../../../functions/containers/csr_matrix/CsrMatrix/nnz.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/tests/cpp_tests/containers/test_csr_matrix/KRATOS_TEST_CASE_IN_SUITE.md)