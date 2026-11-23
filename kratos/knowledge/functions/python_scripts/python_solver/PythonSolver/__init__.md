---
type: Python Method
title: __init__
resource: python_scripts/python_solver.py#L10-L33
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/python_solver/PythonSolver/ValidateSettings
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __init__(self, model, settings):`

# Calls

- [ValidateSettings](../../../../functions/python_scripts/python_solver/PythonSolver/ValidateSettings.md)
- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)