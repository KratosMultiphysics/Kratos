---
type: C++ Function
title: ModelPartOperation
resource: utilities/model_part_operation_utilities.h#L199-L248
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/FillNodesPointerSet
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/FillNodePointersForEntities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/AddNeighbours
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void ModelPartOperation( std::vector<ModelPart::NodeType*>& rOutputNodes, std::vector<ModelPart::ConditionType*>& rOutputConditions, std::vector<ModelPart::ElementType*>& rOutputElements, ModelPart& rMainModelPart, const std::vector<ModelPart const*>& rOperands, const bool AddNeighbourEntities)`

# Calls

- [FillNodesPointerSet](../../../functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/FillNodesPointerSet.md)
- [FillNodePointersForEntities](../../../functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/FillNodePointersForEntities.md)
- [AddNeighbours](../../../functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/AddNeighbours.md)