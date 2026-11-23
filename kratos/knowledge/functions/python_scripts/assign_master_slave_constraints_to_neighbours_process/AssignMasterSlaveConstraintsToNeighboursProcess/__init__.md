---
type: Python Method
title: __init__
resource: python_scripts/assign_master_slave_constraints_to_neighbours_process.py#L21-L101
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetStringArray
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __init__(self, model, settings):`

# Calls

- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [GetDouble](../../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)
- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [GetStringArray](../../../../functions/sources/kratos_parameters/Parameters/GetStringArray.md)
- [GetVariable](../../../../functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable.md)