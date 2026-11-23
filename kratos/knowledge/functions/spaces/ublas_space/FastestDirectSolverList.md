---
type: C++ Function
title: FastestDirectSolverList
resource: spaces/ublas_space.h#L927-L936
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/ConstructLinearSolverFromSettings
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/python_linear_solver_factory/CreateFastestAvailableDirectLinearSolver
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`inline static std::vector<std::string> FastestDirectSolverList()`

# Called by

- [ConstructLinearSolverFromSettings](../../../functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/ConstructLinearSolverFromSettings.md)
- [CreateFastestAvailableDirectLinearSolver](../../../functions/python_scripts/python_linear_solver_factory/CreateFastestAvailableDirectLinearSolver.md)