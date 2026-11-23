---
type: C++ Method
title: ConstructLinearSolverFromSettings
resource: linear_solvers/fallback_linear_solver.h#L630-L648
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spaces/ublas_space/FastestDirectSolverList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetString
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/FallbackLinearSolver
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/AddSolver
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`LinearSolverPointer ConstructLinearSolverFromSettings(const Parameters Settings)`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [FastestDirectSolverList](../../../../functions/spaces/ublas_space/FastestDirectSolverList.md)
- [SetString](../../../../functions/sources/kratos_parameters/Parameters/SetString.md)

# Called by

- [FallbackLinearSolver](../../../../functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/FallbackLinearSolver.md)
- [AddSolver](../../../../functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/AddSolver.md)