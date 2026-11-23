---
type: C++ Method
title: SetSolvers
resource: linear_solvers/fallback_linear_solver.h#L772-L785
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/RemoveValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddEmptyArray
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/FillParametersFromSolver
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SetSolvers(const std::vector<LinearSolverPointer>& rSolvers)`

# Calls

- [RemoveValue](../../../../functions/sources/kratos_parameters/Parameters/RemoveValue.md)
- [AddEmptyArray](../../../../functions/sources/kratos_parameters/Parameters/AddEmptyArray.md)
- [FillParametersFromSolver](../../../../functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/FillParametersFromSolver.md)