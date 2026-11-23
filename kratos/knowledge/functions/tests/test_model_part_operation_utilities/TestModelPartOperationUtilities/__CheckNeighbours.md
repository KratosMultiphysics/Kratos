---
type: Python Method
title: __CheckNeighbours
resource: tests/test_model_part_operation_utilities.py#L36-L57
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/test_model_part_operation_utilities/TestModelPartOperationUtilities/__GenerateUniqueIdsList
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_model_part_operation_utilities/TestModelPartOperationUtilities/test_Union
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_model_part_operation_utilities/TestModelPartOperationUtilities/test_Substract
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_model_part_operation_utilities/TestModelPartOperationUtilities/test_Intersect
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __CheckNeighbours(self, node_ids: 'list[int]', element_ids: 'list[int]', merged_model_part: Kratos.ModelPart, merged_model_part_with_neighbours: Kratos.ModelPart):`

# Calls

- [__GenerateUniqueIdsList](../../../../functions/tests/test_model_part_operation_utilities/TestModelPartOperationUtilities/__GenerateUniqueIdsList.md)

# Called by

- [test_Union](../../../../functions/tests/test_model_part_operation_utilities/TestModelPartOperationUtilities/test_Union.md)
- [test_Substract](../../../../functions/tests/test_model_part_operation_utilities/TestModelPartOperationUtilities/test_Substract.md)
- [test_Intersect](../../../../functions/tests/test_model_part_operation_utilities/TestModelPartOperationUtilities/test_Intersect.md)