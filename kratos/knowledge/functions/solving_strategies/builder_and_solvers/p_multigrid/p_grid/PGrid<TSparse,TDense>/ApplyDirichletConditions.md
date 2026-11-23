---
type: C++ Method
title: ApplyDirichletConditions
resource: solving_strategies/builder_and_solvers/p_multigrid/p_grid.cpp#L301-L358
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/sparse_matrix_multiplication_utility/TransposeMatrix
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void PGrid<TSparse,TDense>::ApplyDirichletConditions(typename IndirectDofSet::const_iterator itParentDofBegin, [[maybe_unused]] typename IndirectDofSet::const_iterator itParentDofEnd)`

# Calls

- [TransposeMatrix](../../../../../../functions/utilities/sparse_matrix_multiplication_utility/TransposeMatrix.md)