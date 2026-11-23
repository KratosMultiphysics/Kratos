---
type: Python Function
title: ConstructSolver
resource: python_scripts/eigen_solver_factory.py#L7-L77
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/kratos_utilities/CheckIfApplicationsAvailable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/python_linear_solver_factory/CreateFastestAvailableDirectLinearSolver
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def ConstructSolver(settings):`

# Calls

- [GetString](../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [CheckIfApplicationsAvailable](../../../functions/python_scripts/kratos_utilities/CheckIfApplicationsAvailable.md)
- [GetBool](../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [CreateFastestAvailableDirectLinearSolver](../../../functions/python_scripts/python_linear_solver_factory/CreateFastestAvailableDirectLinearSolver.md)