---
type: Python Method
title: __init__
resource: python_scripts/sub_model_part_entities_boolean_operation_process.py#L23-L40
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
  - target: functions/python_scripts/sub_model_part_entities_boolean_operation_process/SubModelPartEntitiesBooleanOperationProcess/_GetOrCreateSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/sub_model_part_entities_boolean_operation_process/SubModelPartEntitiesBooleanOperationProcess/_CheckInput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/sub_model_part_entities_boolean_operation_process/SubModelPartEntitiesBooleanOperationProcess/_ExecuteBooleanOperation
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __init__(self, model, settings):`

# Calls

- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [_GetOrCreateSubModelPart](../../../../functions/python_scripts/sub_model_part_entities_boolean_operation_process/SubModelPartEntitiesBooleanOperationProcess/_GetOrCreateSubModelPart.md)
- [_CheckInput](../../../../functions/python_scripts/sub_model_part_entities_boolean_operation_process/SubModelPartEntitiesBooleanOperationProcess/_CheckInput.md)
- [_ExecuteBooleanOperation](../../../../functions/python_scripts/sub_model_part_entities_boolean_operation_process/SubModelPartEntitiesBooleanOperationProcess/_ExecuteBooleanOperation.md)