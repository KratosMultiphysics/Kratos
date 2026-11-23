---
type: C++ Function
title: Union
resource: utilities/sub_model_part_entities_boolean_operation_utility.h#L76-L79
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
  - target: functions/tests/test_model_part_operation_utilities/TestModelPartOperationUtilities/test_Order
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_model_part_operation_utilities/TestModelPartOperationUtilities/test_Union
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_model_part_operation_utilities/TestModelPartOperationUtilities/test_Sum
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void Union(ModelPart& rModelPartA, ModelPart& rModelPartB, ModelPart& rDestination)`

# Calls

- [BooleanOperation](../../../functions/utilities/sub_model_part_entities_boolean_operation_utility/SubModelPartEntitiesBooleanOperationUtility<
    TEntityType, TContainerType>/BooleanOperation.md)

# Called by

- [_ExecuteBooleanOperation](../../../functions/python_scripts/sub_model_part_entities_boolean_operation_process/SubModelPartEntitiesBooleanOperationProcess/_ExecuteBooleanOperation.md)
- [test_Order](../../../functions/tests/test_model_part_operation_utilities/TestModelPartOperationUtilities/test_Order.md)
- [test_Union](../../../functions/tests/test_model_part_operation_utilities/TestModelPartOperationUtilities/test_Union.md)
- [test_Sum](../../../functions/tests/test_model_part_operation_utilities/TestModelPartOperationUtilities/test_Sum.md)