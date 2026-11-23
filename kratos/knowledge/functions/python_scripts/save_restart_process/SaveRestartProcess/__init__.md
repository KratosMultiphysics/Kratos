---
type: Python Method
title: __init__
resource: python_scripts/save_restart_process.py#L15-L59
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/AddValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/RemoveValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/restart_utility/RestartUtility/CreateOutputFolder
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __init__(self, model, params):`

# Calls

- [AddValue](../../../../functions/sources/kratos_parameters/Parameters/AddValue.md)
- [RemoveValue](../../../../functions/sources/kratos_parameters/Parameters/RemoveValue.md)
- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [CreateOutputFolder](../../../../functions/python_scripts/restart_utility/RestartUtility/CreateOutputFolder.md)