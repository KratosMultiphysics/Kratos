---
type: C++ Method
title: FillParametersFromSolver
resource: linear_solvers/fallback_linear_solver.h#L657-L675
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
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
  - target: functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/SetSolvers
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void FillParametersFromSolver(LinearSolverPointer pSolver)`

# Calls

- [SetString](../../../../functions/sources/kratos_parameters/Parameters/SetString.md)

# Called by

- [FallbackLinearSolver](../../../../functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/FallbackLinearSolver.md)
- [AddSolver](../../../../functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/AddSolver.md)
- [SetSolvers](../../../../functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/SetSolvers.md)