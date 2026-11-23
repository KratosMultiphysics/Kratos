---
type: Python Method
title: __init__
resource: python_scripts/assign_time_derivative_process.py#L16-L33
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/RemoveValue
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __init__(self, Model, settings ):`

# Calls

- [GetVariable](../../../../functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [RemoveValue](../../../../functions/sources/kratos_parameters/Parameters/RemoveValue.md)