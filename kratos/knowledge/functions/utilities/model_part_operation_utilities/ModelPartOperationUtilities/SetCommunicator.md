---
type: C++ Method
title: SetCommunicator
resource: utilities/model_part_operation_utilities.cpp#L151-L178
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/communicator/Communicator/SetNumberOfColors
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/GetNumberOfColors
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/FillWithMainNodesFromSearchNodes
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartOperationUtilities::SetCommunicator( ModelPart& rOutputModelPart, ModelPart& rMainModelPart)`

# Calls

- [SetNumberOfColors](../../../../functions/sources/communicator/Communicator/SetNumberOfColors.md)
- [GetNumberOfColors](../../../../functions/sources/communicator/Communicator/GetNumberOfColors.md)
- [FillWithMainNodesFromSearchNodes](../../../../functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/FillWithMainNodesFromSearchNodes.md)