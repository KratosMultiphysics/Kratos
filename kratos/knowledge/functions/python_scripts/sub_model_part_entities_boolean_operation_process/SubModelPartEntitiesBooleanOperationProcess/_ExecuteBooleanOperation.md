---
type: Python Method
title: _ExecuteBooleanOperation
resource: python_scripts/sub_model_part_entities_boolean_operation_process.py#L61-L75
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/sub_model_part_entities_boolean_operation_utility/Union
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/sub_model_part_entities_boolean_operation_utility/Difference
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/sub_model_part_entities_boolean_operation_process/SubModelPartEntitiesBooleanOperationProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _ExecuteBooleanOperation(self):`

# Calls

- [Union](../../../../functions/utilities/sub_model_part_entities_boolean_operation_utility/Union.md)
- [Difference](../../../../functions/utilities/sub_model_part_entities_boolean_operation_utility/Difference.md)

# Called by

- [__init__](../../../../functions/python_scripts/sub_model_part_entities_boolean_operation_process/SubModelPartEntitiesBooleanOperationProcess/__init__.md)