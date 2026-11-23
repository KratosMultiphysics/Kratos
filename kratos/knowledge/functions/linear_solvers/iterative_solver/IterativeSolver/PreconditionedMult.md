---
type: C++ Method
title: PreconditionedMult
resource: linear_solvers/iterative_solver.h#L416-L419
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/linear_solvers/bicgstab_solver/BICGSTABSolver/IterativeSolve
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/linear_solvers/cg_solver/CGSolver/IterativeSolve
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/linear_solvers/deflated_cg_solver/DeflatedCGSolver/IterativeSolve
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/linear_solvers/tfqmr_solver/TFQMRSolver/IterativeSolve
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void PreconditionedMult(SparseMatrixType& rA, VectorType& rX, VectorType& rY)`

# Called by

- [IterativeSolve](../../../../functions/linear_solvers/bicgstab_solver/BICGSTABSolver/IterativeSolve.md)
- [IterativeSolve](../../../../functions/linear_solvers/cg_solver/CGSolver/IterativeSolve.md)
- [IterativeSolve](../../../../functions/linear_solvers/deflated_cg_solver/DeflatedCGSolver/IterativeSolve.md)
- [IterativeSolve](../../../../functions/linear_solvers/tfqmr_solver/TFQMRSolver/IterativeSolve.md)