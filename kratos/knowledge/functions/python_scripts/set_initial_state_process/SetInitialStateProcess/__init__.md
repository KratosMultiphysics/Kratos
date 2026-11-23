---
type: Python Method
title: __init__
resource: python_scripts/set_initial_state_process.py#L57-L107
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/interval_utility/IntervalUtility/IntervalUtility
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/set_initial_state_process/components_to_functions
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __init__(self, Model, settings):`

# Calls

- [IntervalUtility](../../../../functions/utilities/interval_utility/IntervalUtility/IntervalUtility.md)
- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [components_to_functions](../../../../functions/python_scripts/set_initial_state_process/components_to_functions.md)