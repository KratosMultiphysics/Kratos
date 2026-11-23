---
type: C++ Function
title: InPlaceMatrixAdd
resource: solving_strategies/builder_and_solvers/p_multigrid/sparse_utilities.hpp#L79-L118
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/csr_matrix/CsrMatrix/nnz
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler/AugmentedLagrangeConstraintAssembler<TSparse,TDense>/Initialize
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void InPlaceMatrixAdd(typename TUblasSparseSpace<TValue>::MatrixType& rLeft, const typename TUblasSparseSpace<TValue>::MatrixType& rRight, const TValue Coefficient = 1.0)`

# Calls

- [nnz](../../../../../functions/containers/csr_matrix/CsrMatrix/nnz.md)

# Called by

- [Initialize](../../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler/AugmentedLagrangeConstraintAssembler<TSparse,TDense>/Initialize.md)