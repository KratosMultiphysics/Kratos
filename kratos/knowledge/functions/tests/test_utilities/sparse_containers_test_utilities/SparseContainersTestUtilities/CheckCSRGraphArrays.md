---
type: C++ Method
title: CheckCSRGraphArrays
resource: tests/test_utilities/sparse_containers_test_utilities.cpp#L158-L179
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/tests/cpp_tests/containers/test_sparse_graph/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool SparseContainersTestUtilities::CheckCSRGraphArrays( const Kratos::span<IndexType>& rRowIndices, const Kratos::span<IndexType>& rColIndices, const MatrixMapType& rReferenceGraph)`

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/tests/cpp_tests/containers/test_sparse_graph/KRATOS_TEST_CASE_IN_SUITE.md)