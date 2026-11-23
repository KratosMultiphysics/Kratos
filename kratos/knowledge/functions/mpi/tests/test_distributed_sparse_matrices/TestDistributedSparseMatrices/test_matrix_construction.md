---
type: Python Method
title: test_matrix_construction
resource: mpi/tests/test_distributed_sparse_matrices.py#L83-L190
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/mpi/tests/test_distributed_sparse_matrices/TestDistributedSparseMatrices/GetLocalConnectivities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/distributed_numbering/DistributedNumbering/GlobalId
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/Matrix/fill
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/distributed_vector_importer/DistributedVectorImporter/ImportData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/dot
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/distributed_csr_matrix/DistributedCsrMatrix/ToSerialCSR
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/KratosUnittest/TestCase/assertVectorAlmostEqual
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def test_matrix_construction(self):`

# Calls

- [GetLocalConnectivities](../../../../../functions/mpi/tests/test_distributed_sparse_matrices/TestDistributedSparseMatrices/GetLocalConnectivities.md)
- [GlobalId](../../../../../functions/containers/distributed_numbering/DistributedNumbering/GlobalId.md)
- [fill](../../../../../functions/includes/amatrix_interface/Matrix/fill.md)
- [ImportData](../../../../../functions/containers/distributed_vector_importer/DistributedVectorImporter/ImportData.md)
- [dot](../../../../../functions/includes/amatrix_interface/dot.md)
- [ToSerialCSR](../../../../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/ToSerialCSR.md)
- [assertVectorAlmostEqual](../../../../../functions/python_scripts/KratosUnittest/TestCase/assertVectorAlmostEqual.md)