---
type: C++ Function
title: FillDeflatedMatrix
resource: utilities/deflation_utils.h#L290-L490
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/linear_solvers/deflated_cg_solver/DeflatedCGSolver/IterativeSolve
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void FillDeflatedMatrix( const SparseMatrixType& rA, std::vector<int>& w, SparseMatrixType& Ah)`

# Called by

- [IterativeSolve](../../../functions/linear_solvers/deflated_cg_solver/DeflatedCGSolver/IterativeSolve.md)