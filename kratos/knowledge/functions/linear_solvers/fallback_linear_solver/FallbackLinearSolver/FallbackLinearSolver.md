---
type: C++ Method
title: FallbackLinearSolver
resource: linear_solvers/fallback_linear_solver.h#L110
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/ConstructLinearSolverFromSettings
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/CommonSettingsFromParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/FillParametersFromSolver
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`explicit FallbackLinearSolver() = default;`

# Calls

- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [ConstructLinearSolverFromSettings](../../../../functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/ConstructLinearSolverFromSettings.md)
- [CommonSettingsFromParameters](../../../../functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/CommonSettingsFromParameters.md)
- [FillParametersFromSolver](../../../../functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/FillParametersFromSolver.md)