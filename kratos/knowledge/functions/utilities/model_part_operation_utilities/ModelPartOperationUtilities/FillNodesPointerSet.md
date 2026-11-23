---
type: C++ Method
title: FillNodesPointerSet
resource: utilities/model_part_operation_utilities.cpp#L31-L38
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/FillWithMainNodesFromSearchNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/CheckValidityOfModelPartsForOperations
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/model_part_operation_utilities/ModelPartOperation
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartOperationUtilities::FillNodesPointerSet( std::set<ModelPart::NodeType const*>& rOutput, const ModelPart::NodesContainerType& rNodes)`

# Called by

- [FillWithMainNodesFromSearchNodes](../../../../functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/FillWithMainNodesFromSearchNodes.md)
- [CheckValidityOfModelPartsForOperations](../../../../functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/CheckValidityOfModelPartsForOperations.md)
- [ModelPartOperation](../../../../functions/utilities/model_part_operation_utilities/ModelPartOperation.md)