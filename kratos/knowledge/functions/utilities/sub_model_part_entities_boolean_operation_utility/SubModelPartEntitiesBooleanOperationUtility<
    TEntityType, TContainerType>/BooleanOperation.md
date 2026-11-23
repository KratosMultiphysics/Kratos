---
type: C++ Method
title: BooleanOperation
resource: utilities/sub_model_part_entities_boolean_operation_utility.cpp#L101-L139
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/IsSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: |-
      functions/utilities/sub_model_part_entities_boolean_operation_utility/SubModelPartEntitiesBooleanOperationUtility<
          TEntityType, TContainerType>/GetContainerIds
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/sub_model_part_entities_boolean_operation_utility/Union
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/sub_model_part_entities_boolean_operation_utility/Intersection
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/sub_model_part_entities_boolean_operation_utility/Difference
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SubModelPartEntitiesBooleanOperationUtility< TEntityType, TContainerType>::BooleanOperation( ModelPart& rModelPartA, ModelPart& rModelPartB, ModelPart& rDestination, BooleanOperators ThisOperator)`

# Calls

- [IsSubModelPart](../../../../functions/includes/model_part/IsSubModelPart.md)
- [GetContainerIds](../../../../functions/utilities/sub_model_part_entities_boolean_operation_utility/SubModelPartEntitiesBooleanOperationUtility<
    TEntityType, TContainerType>/GetContainerIds.md)

# Called by

- [Union](../../../../functions/utilities/sub_model_part_entities_boolean_operation_utility/Union.md)
- [Intersection](../../../../functions/utilities/sub_model_part_entities_boolean_operation_utility/Intersection.md)
- [Difference](../../../../functions/utilities/sub_model_part_entities_boolean_operation_utility/Difference.md)