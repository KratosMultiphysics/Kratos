---
type: Python Method
title: __init__
resource: python_scripts/assign_scalar_input_to_entities_process.py#L20-L86
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/interval_utility/IntervalUtility/IntervalUtility
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/IsString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __init__(self, Model, settings ):`

# Calls

- [IntervalUtility](../../../../functions/utilities/interval_utility/IntervalUtility/IntervalUtility.md)
- [IsString](../../../../functions/sources/kratos_parameters/Parameters/IsString.md)
- [SetString](../../../../functions/sources/kratos_parameters/Parameters/SetString.md)
- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [AddValue](../../../../functions/sources/kratos_parameters/Parameters/AddValue.md)
- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)