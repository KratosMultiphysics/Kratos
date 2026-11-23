---
type: Python Method
title: ValidateSettings
resource: python_scripts/python_solver.py#L46-L50
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/python_solver/PythonSolver/__init__
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def ValidateSettings(self):`

# Calls

- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)

# Called by

- [__init__](../../../../functions/python_scripts/python_solver/PythonSolver/__init__.md)