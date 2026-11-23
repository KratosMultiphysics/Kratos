---
type: C++ Module
title: test_distributed_csr_matrix
resource: mpi/tests/cpp_tests/containers/test_distributed_csr_matrix.cpp#L1-L468
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/containers-sparse-contiguous-row-graph-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/containers-distributed-sparse-graph-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/containers-distributed-csr-matrix-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/containers-distributed-system-vector-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/containers-distributed-vector-importer-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/mpi-includes-mpi-data-communicator-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/mpi-testing-mpi-testing-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/mpi-tests-test-utilities-distributed-sparse-containers-test-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/mpi-utilities-amgcl-distributed-csr-conversion-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/mpi-utilities-amgcl-distributed-csr-spmm-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-openmp-utils-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-builtin-timer-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/mpi/tests/cpp_tests/containers/test_distributed_csr_matrix/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/mpi/tests/cpp_tests/containers/test_distributed_csr_matrix/KRATOS_TEST_CASE_IN_SUITE-2.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/mpi/tests/cpp_tests/containers/test_distributed_csr_matrix/KRATOS_TEST_CASE_IN_SUITE-3.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/mpi/tests/cpp_tests/containers/test_distributed_csr_matrix/KRATOS_TEST_CASE_IN_SUITE-4.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/mpi/tests/cpp_tests/containers/test_distributed_csr_matrix/KRATOS_TEST_CASE_IN_SUITE-5.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/mpi/tests/cpp_tests/containers/test_distributed_csr_matrix/KRATOS_TEST_CASE_IN_SUITE-6.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/mpi/tests/cpp_tests/containers/test_distributed_csr_matrix/KRATOS_TEST_CASE_IN_SUITE-7.md)

# Imports

- `"containers/sparse_contiguous_row_graph.h"`
- `"containers/distributed_sparse_graph.h"`
- `"containers/distributed_csr_matrix.h"`
- `"containers/distributed_system_vector.h"`
- `"containers/distributed_vector_importer.h"`
- `"mpi/includes/mpi_data_communicator.h"`
- `"mpi/testing/mpi_testing.h"`
- `"mpi/tests/test_utilities/distributed_sparse_containers_test_utilities.h"`
- `"mpi/utilities/amgcl_distributed_csr_conversion_utilities.h"`
- `"mpi/utilities/amgcl_distributed_csr_spmm_utilities.h"`
- `"utilities/openmp_utils.h"`
- `"utilities/builtin_timer.h"`

# Member of

- [KratosMultiphysics](../../../../../packages/KratosMultiphysics.md)