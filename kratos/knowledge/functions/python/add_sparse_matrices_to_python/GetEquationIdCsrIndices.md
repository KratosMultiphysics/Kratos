---
type: C++ Function
title: GetEquationIdCsrIndices
resource: python/add_sparse_matrices_to_python.cpp#L74-L105
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/csr_matrix/CsrMatrix/FindValueIndex
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`pybind11::array_t<TIndexType> GetEquationIdCsrIndices( const CsrMatrixType& rCsrMatrix, const pybind11::array_t<TIndexType>& rConnectivities)`

# Calls

- [FindValueIndex](../../../functions/containers/csr_matrix/CsrMatrix/FindValueIndex.md)