---
type: Python Function
title: Factory
resource: python_scripts/csv_points_output_process.py#L8-L57
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
  - target: functions/sources/kratos_parameters/Parameters/RemoveValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddEmptyValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetMatrix
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def Factory(settings, Model):`

# Calls

- [ValidateAndAssignDefaults](../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [GetString](../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [RemoveValue](../../../functions/sources/kratos_parameters/Parameters/RemoveValue.md)
- [AddEmptyValue](../../../functions/sources/kratos_parameters/Parameters/AddEmptyValue.md)
- [SetMatrix](../../../functions/sources/kratos_parameters/Parameters/SetMatrix.md)