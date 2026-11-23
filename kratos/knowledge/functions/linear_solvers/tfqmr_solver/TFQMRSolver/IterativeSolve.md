---
type: C++ Method
title: IterativeSolve
resource: linear_solvers/tfqmr_solver.h#L191-L361
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/linear_solvers/iterative_solver/IterativeSolver/PreconditionedMult
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool IterativeSolve(SparseMatrixType& rA, VectorType& rX, VectorType& rB)`

# Calls

- [PreconditionedMult](../../../../functions/linear_solvers/iterative_solver/IterativeSolver/PreconditionedMult.md)