---
type: C++ Method
title: DuplicateCommunicatorData
resource: modeler/connectivity_preserve_modeler.cpp#L299-L342
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
  - target: functions/sources/communicator/Communicator/SetLocalMesh
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ConnectivityPreserveModeler::DuplicateCommunicatorData( ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart) const`

# Calls

- [SetNumberOfColors](../../../../functions/sources/communicator/Communicator/SetNumberOfColors.md)
- [GetNumberOfColors](../../../../functions/sources/communicator/Communicator/GetNumberOfColors.md)
- [SetLocalMesh](../../../../functions/sources/communicator/Communicator/SetLocalMesh.md)