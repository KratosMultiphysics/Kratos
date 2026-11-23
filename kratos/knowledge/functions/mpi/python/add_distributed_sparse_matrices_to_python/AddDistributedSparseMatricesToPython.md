---
type: C++ Function
title: AddDistributedSparseMatricesToPython
resource: mpi/python/add_distributed_sparse_matrices_to_python.cpp#L30-L245
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/distributed_csr_matrix/DistributedCsrMatrix/GetDiagonalIndex2DataInGlobalNumbering
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/distributed_csr_matrix/DistributedCsrMatrix/GetOffDiagonalIndex2DataInGlobalNumbering
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/python/kratos_mpi_python/PYBIND11_MODULE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AddDistributedSparseMatricesToPython(pybind11::module& m)`

# Calls

- [GetDiagonalIndex2DataInGlobalNumbering](../../../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/GetDiagonalIndex2DataInGlobalNumbering.md)
- [GetOffDiagonalIndex2DataInGlobalNumbering](../../../../functions/containers/distributed_csr_matrix/DistributedCsrMatrix/GetOffDiagonalIndex2DataInGlobalNumbering.md)

# Called by

- [PYBIND11_MODULE](../../../../functions/mpi/python/kratos_mpi_python/PYBIND11_MODULE.md)