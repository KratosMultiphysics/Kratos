---
type: C++ Function
title: ApplyWtranspose
resource: utilities/deflation_utils.h#L275-L286
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/linear_solvers/deflated_cg_solver/DeflatedCGSolver/IterativeSolve
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void ApplyWtranspose(const std::vector<int>& w, const SparseVectorType& x, SparseVectorType& y)`

# Called by

- [IterativeSolve](../../../functions/linear_solvers/deflated_cg_solver/DeflatedCGSolver/IterativeSolve.md)