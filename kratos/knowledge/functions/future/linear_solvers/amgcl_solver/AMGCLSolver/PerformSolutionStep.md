---
type: C++ Method
title: PerformSolutionStep
resource: future/linear_solvers/amgcl_solver.h#L332-L443
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/future/containers/linear_system/LinearSystem/pGetVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/linear_solvers/amgcl_solver_impl/AMGCLSolve
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool PerformSolutionStep(LinearSystemType& rLinearSystem) override`

# Calls

- [pGetVector](../../../../../functions/future/containers/linear_system/LinearSystem/pGetVector.md)
- [AMGCLSolve](../../../../../functions/future/linear_solvers/amgcl_solver_impl/AMGCLSolve.md)