---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: mpi/tests/cpp_tests/containers/test_distributed_system_vector.cpp#L30-L74
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/mpi/tests/test_utilities/distributed_sparse_containers_test_utilities/DistributedSparseContainersTestUtilities/GetReferenceVectorAsMap
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/distributed_numbering/DistributedNumbering/GlobalId
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/distributed_vector_importer/DistributedVectorImporter/ImportData
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(DistributedSystemVectorAssembly, KratosMPICoreFastSuite)`

# Calls

- [GetReferenceVectorAsMap](../../../../../../functions/mpi/tests/test_utilities/distributed_sparse_containers_test_utilities/DistributedSparseContainersTestUtilities/GetReferenceVectorAsMap.md)
- [GlobalId](../../../../../../functions/containers/distributed_numbering/DistributedNumbering/GlobalId.md)
- [ImportData](../../../../../../functions/containers/distributed_vector_importer/DistributedVectorImporter/ImportData.md)