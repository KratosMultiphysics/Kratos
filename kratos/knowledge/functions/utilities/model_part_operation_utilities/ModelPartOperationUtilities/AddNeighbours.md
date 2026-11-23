---
type: C++ Method
title: AddNeighbours
resource: utilities/model_part_operation_utilities.cpp#L122-L139
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/variable_utils/SetFlag
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/FindNeighbourEntities
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/model_part_operation_utilities/ModelPartOperation
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartOperationUtilities::AddNeighbours( std::vector<ModelPart::NodeType*>& rOutputNodes, std::vector<ModelPart::ConditionType*>& rOutputConditions, std::vector<ModelPart::ElementType*>& rOutputElements, ModelPart& rMainModelPart)`

# Calls

- [SetFlag](../../../../functions/utilities/variable_utils/SetFlag.md)
- [FindNeighbourEntities](../../../../functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/FindNeighbourEntities.md)

# Called by

- [ModelPartOperation](../../../../functions/utilities/model_part_operation_utilities/ModelPartOperation.md)