---
type: C++ Method
title: SetPreconditioner
resource: linear_solvers/iterative_solver.h#L270-L273
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/linear_solvers/bicgstab_solver/BICGSTABSolver/BICGSTABSolver
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/linear_solvers/cg_solver/CGSolver/CGSolver
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/linear_solvers/tfqmr_solver/TFQMRSolver/TFQMRSolver
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual void SetPreconditioner(typename TPreconditionerType::Pointer pNewPreconditioner)`

# Called by

- [BICGSTABSolver](../../../../functions/linear_solvers/bicgstab_solver/BICGSTABSolver/BICGSTABSolver.md)
- [CGSolver](../../../../functions/linear_solvers/cg_solver/CGSolver/CGSolver.md)
- [TFQMRSolver](../../../../functions/linear_solvers/tfqmr_solver/TFQMRSolver/TFQMRSolver.md)