---
type: C++ Method
title: GetOffDiagonalIndex2DataInGlobalNumbering
resource: containers/distributed_csr_matrix.h#L431-L439
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/distributed_csr_matrix/DistributedCsrMatrix/GetOffDiaGlobalId
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/containers/distributed_csr_matrix/DistributedCsrMatrix/PrintData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/python/add_distributed_sparse_matrices_to_python/AddDistributedSparseMatricesToPython
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/cpp_tests/containers/test_distributed_csr_matrix/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/amgcl_distributed_csr_conversion_utilities/Transpose
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/amgcl_distributed_csr_spmm_utilities/SparseMultiply
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/petsc_conversion_tools/PetscCSRWrapper/__init__
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`DenseVector<TIndexType> GetOffDiagonalIndex2DataInGlobalNumbering() const`

# Calls

- [GetOffDiaGlobalId](../../../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/GetOffDiaGlobalId.md)

# Called by

- [PrintData](../../../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/PrintData.md)
- [AddDistributedSparseMatricesToPython](../../../../functions/mpi/python/add_distributed_sparse_matrices_to_python/AddDistributedSparseMatricesToPython.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/mpi/tests/cpp_tests/containers/test_distributed_csr_matrix/KRATOS_TEST_CASE_IN_SUITE.md)
- [Transpose](../../../../functions/mpi/utilities/amgcl_distributed_csr_conversion_utilities/Transpose.md)
- [SparseMultiply](../../../../functions/mpi/utilities/amgcl_distributed_csr_spmm_utilities/SparseMultiply.md)
- [__init__](../../../../functions/python_scripts/petsc_conversion_tools/PetscCSRWrapper/__init__.md)