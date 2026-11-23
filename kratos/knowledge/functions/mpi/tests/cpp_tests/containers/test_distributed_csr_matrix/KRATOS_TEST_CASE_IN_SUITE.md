---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: mpi/tests/cpp_tests/containers/test_distributed_csr_matrix.cpp#L33-L68
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/distributed_numbering/DistributedNumbering/GlobalId
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/distributed_csr_matrix/DistributedCsrMatrix/GetOffDiagonalIndex2DataInGlobalNumbering
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/distributed_numbering/DistributedNumbering/IsLocal
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(DistributedCSRMatrixConstruction, KratosMPICoreFastSuite)`

# Calls

- [GlobalId](../../../../../../functions/containers/distributed_numbering/DistributedNumbering/GlobalId.md)
- [GetOffDiagonalIndex2DataInGlobalNumbering](../../../../../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/GetOffDiagonalIndex2DataInGlobalNumbering.md)
- [IsLocal](../../../../../../functions/containers/distributed_numbering/DistributedNumbering/IsLocal.md)