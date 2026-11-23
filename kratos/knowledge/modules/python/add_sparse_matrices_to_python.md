---
type: C++ Module
title: add_sparse_matrices_to_python
resource: python/add_sparse_matrices_to_python.cpp#L1-L513
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/pybind11-numpy-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-define-python-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-ublas-interface-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/python-add-sparse-matrices-to-python-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/python-numpy-utils-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/containers-nd-data-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/containers-sparse-graph-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/containers-sparse-contiguous-row-graph-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/containers-csr-matrix-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/containers-system-vector-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-amgcl-csr-spmm-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [GetEquationIdCsrIndices](../../functions/python/add_sparse_matrices_to_python/GetEquationIdCsrIndices.md)
- [GetEquationIdCsrIndices](../../functions/python/add_sparse_matrices_to_python/GetEquationIdCsrIndices-2.md)
- [AssembleWithCsrIndicesTemplated](../../functions/python/add_sparse_matrices_to_python/AssembleWithCsrIndicesTemplated.md)
- [AssembleWithEquationIdsTemplated](../../functions/python/add_sparse_matrices_to_python/AssembleWithEquationIdsTemplated.md)
- [AddSparseMatricesToPython](../../functions/python/add_sparse_matrices_to_python/AddSparseMatricesToPython.md)

# Imports

- `<pybind11/numpy.h>`
- `"includes/define_python.h"`
- `"includes/ublas_interface.h"`
- `"python/add_sparse_matrices_to_python.h"`
- `"python/numpy_utils.h"`
- `"containers/nd_data.h"`
- `"containers/sparse_graph.h"`
- `"containers/sparse_contiguous_row_graph.h"`
- `"containers/csr_matrix.h"`
- `"containers/system_vector.h"`
- `"utilities/amgcl_csr_spmm_utilities.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)