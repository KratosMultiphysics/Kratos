---
type: C++ Function
title: MergeMatrices
resource: solving_strategies/builder_and_solvers/p_multigrid/sparse_utilities.hpp#L139-L229
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/csr_matrix/CsrMatrix/nnz
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/pointer_vector/emplace_back
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void MergeMatrices(typename TUblasSparseSpace<TValue>::MatrixType& rLeft, const typename TUblasSparseSpace<TValue>::MatrixType& rRight)`

# Calls

- [nnz](../../../../../functions/containers/csr_matrix/CsrMatrix/nnz.md)
- [emplace_back](../../../../../functions/containers/pointer_vector/emplace_back.md)