---
type: C++ Method
title: CGSolver
resource: future/linear_solvers/cg_solver.h#L82-L97
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/linear_solvers/iterative_solver/IterativeSolver/SetLeftPreconditioner
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`CGSolver(Parameters Settings = Parameters(R"({})")) : BaseType()`

# Calls

- [ValidateAndAssignDefaults](../../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [GetString](../../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [SetLeftPreconditioner](../../../../../functions/future/linear_solvers/iterative_solver/IterativeSolver/SetLeftPreconditioner.md)