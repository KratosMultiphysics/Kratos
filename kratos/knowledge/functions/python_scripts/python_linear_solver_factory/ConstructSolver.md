---
type: Python Function
title: ConstructSolver
resource: python_scripts/python_linear_solver_factory.py#L5-L20
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/python_linear_solver_factory/__GetSolverTypeAndImportApplication
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def ConstructSolver(configuration):`

# Calls

- [__GetSolverTypeAndImportApplication](../../../functions/python_scripts/python_linear_solver_factory/__GetSolverTypeAndImportApplication.md)
- [GetString](../../../functions/sources/kratos_parameters/Parameters/GetString.md)