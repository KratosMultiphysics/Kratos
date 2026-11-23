---
type: Python Method
title: _GetRestartSettings
resource: python_scripts/python_solver.py#L169-L177
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/RemoveValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddValue
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/python_solver/PythonSolver/_ImportModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _GetRestartSettings(self, model_part_import_settings):`

# Calls

- [RemoveValue](../../../../functions/sources/kratos_parameters/Parameters/RemoveValue.md)
- [AddValue](../../../../functions/sources/kratos_parameters/Parameters/AddValue.md)

# Called by

- [_ImportModelPart](../../../../functions/python_scripts/python_solver/PythonSolver/_ImportModelPart.md)