---
type: Python Method
title: __init__
resource: python_scripts/output_process_with_controller.py#L18-L26
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/model_parameters_factory/KratosModelParametersFactory/ConstructItem
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __init__(self, model, params):`

# Calls

- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [ConstructItem](../../../../functions/python_scripts/model_parameters_factory/KratosModelParametersFactory/ConstructItem.md)