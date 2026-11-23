---
type: Python Method
title: __init__
resource: python_scripts/multiple_points_output_process.py#L20-L64
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/RemoveValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddEmptyValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetVector
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __init__(self, model, params):`

# Calls

- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [RemoveValue](../../../../functions/sources/kratos_parameters/Parameters/RemoveValue.md)
- [AddEmptyValue](../../../../functions/sources/kratos_parameters/Parameters/AddEmptyValue.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [SetString](../../../../functions/sources/kratos_parameters/Parameters/SetString.md)
- [SetVector](../../../../functions/sources/kratos_parameters/Parameters/SetVector.md)