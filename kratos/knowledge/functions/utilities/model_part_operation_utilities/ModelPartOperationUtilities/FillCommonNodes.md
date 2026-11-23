---
type: C++ Method
title: FillCommonNodes
resource: utilities/model_part_operation_utilities.cpp#L62-L82
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/FillWithMainNodesFromSearchNodes
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartOperationUtilities::FillCommonNodes( ModelPart::NodesContainerType& rOutput, ModelPart::NodesContainerType& rMainNodes, const std::set<ModelPart::NodeType const*>& rNodesSet)`

# Called by

- [FillWithMainNodesFromSearchNodes](../../../../functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/FillWithMainNodesFromSearchNodes.md)