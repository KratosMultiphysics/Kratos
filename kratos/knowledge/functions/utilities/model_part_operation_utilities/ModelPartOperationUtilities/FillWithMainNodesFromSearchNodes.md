---
type: C++ Method
title: FillWithMainNodesFromSearchNodes
resource: utilities/model_part_operation_utilities.cpp#L141-L149
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/FillNodesPointerSet
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/FillCommonNodes
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/SetCommunicator
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartOperationUtilities::FillWithMainNodesFromSearchNodes( ModelPart::NodesContainerType& rOutput, ModelPart::NodesContainerType& rMainNodes, ModelPart::NodesContainerType& rSearchedNodes)`

# Calls

- [FillNodesPointerSet](../../../../functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/FillNodesPointerSet.md)
- [FillCommonNodes](../../../../functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/FillCommonNodes.md)

# Called by

- [SetCommunicator](../../../../functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/SetCommunicator.md)