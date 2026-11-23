---
type: C++ Method
title: RandomElementConnectivities
resource: tests/test_utilities/sparse_containers_test_utilities.cpp#L31-L72
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/mpi/tests/cpp_tests/containers/test_distributed_sparse_graph/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/containers/test_sparse_graph/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`SparseContainersTestUtilities::ElementConnectivityType SparseContainersTestUtilities::RandomElementConnectivities( const IndexType BlockSize, const IndexType NodesInElem, const IndexType IndexBegin, const IndexType IndexEnd, const IndexType NDof, const IndexType StandardDev)`

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/mpi/tests/cpp_tests/containers/test_distributed_sparse_graph/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/tests/cpp_tests/containers/test_sparse_graph/KRATOS_TEST_CASE_IN_SUITE.md)