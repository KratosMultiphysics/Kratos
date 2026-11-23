---
type: C++ Function
title: Difference
resource: utilities/sub_model_part_entities_boolean_operation_utility.h#L86-L89
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: |-
      functions/utilities/sub_model_part_entities_boolean_operation_utility/SubModelPartEntitiesBooleanOperationUtility<
          TEntityType, TContainerType>/BooleanOperation
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/sub_model_part_entities_boolean_operation_process/SubModelPartEntitiesBooleanOperationProcess/_ExecuteBooleanOperation
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void Difference(ModelPart& rModelPartA, ModelPart& rModelPartB, ModelPart& rDestination)`

# Calls

- [BooleanOperation](../../../functions/utilities/sub_model_part_entities_boolean_operation_utility/SubModelPartEntitiesBooleanOperationUtility<
    TEntityType, TContainerType>/BooleanOperation.md)

# Called by

- [_ExecuteBooleanOperation](../../../functions/python_scripts/sub_model_part_entities_boolean_operation_process/SubModelPartEntitiesBooleanOperationProcess/_ExecuteBooleanOperation.md)