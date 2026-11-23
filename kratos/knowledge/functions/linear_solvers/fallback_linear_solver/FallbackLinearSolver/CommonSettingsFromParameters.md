---
type: C++ Method
title: CommonSettingsFromParameters
resource: linear_solvers/fallback_linear_solver.h#L684-L688
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/FallbackLinearSolver
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CommonSettingsFromParameters()`

# Calls

- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)

# Called by

- [FallbackLinearSolver](../../../../functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/FallbackLinearSolver.md)