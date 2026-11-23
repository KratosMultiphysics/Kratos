---
type: C++ Method
title: UpdateSolverIndex
resource: linear_solvers/fallback_linear_solver.h#L723-L746
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/GetCurrentSolver
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/PrettyPrintJsonString
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/Solve
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void UpdateSolverIndex()`

# Calls

- [GetCurrentSolver](../../../../functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/GetCurrentSolver.md)
- [PrettyPrintJsonString](../../../../functions/sources/kratos_parameters/Parameters/PrettyPrintJsonString.md)

# Called by

- [Solve](../../../../functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/Solve.md)