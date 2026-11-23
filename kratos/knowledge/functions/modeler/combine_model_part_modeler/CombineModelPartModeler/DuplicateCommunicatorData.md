---
type: C++ Method
title: DuplicateCommunicatorData
resource: modeler/combine_model_part_modeler.cpp#L358-L376
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
  - target: functions/modeler/combine_model_part_modeler/CombineModelPartModeler/PopulateLocalMesh
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CombineModelPartModeler::DuplicateCommunicatorData( ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart) const`

# Calls

- [SetNumberOfColors](../../../../functions/sources/communicator/Communicator/SetNumberOfColors.md)
- [GetNumberOfColors](../../../../functions/sources/communicator/Communicator/GetNumberOfColors.md)
- [PopulateLocalMesh](../../../../functions/modeler/combine_model_part_modeler/CombineModelPartModeler/PopulateLocalMesh.md)