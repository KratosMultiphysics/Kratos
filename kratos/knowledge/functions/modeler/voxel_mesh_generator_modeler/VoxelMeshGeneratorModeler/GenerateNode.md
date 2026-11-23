---
type: C++ Method
title: GenerateNode
resource: modeler/voxel_mesh_generator_modeler.cpp#L262-L272
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/node/SetSolutionStepVariablesList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/pGetNodalSolutionStepVariablesList
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`Node::Pointer VoxelMeshGeneratorModeler::GenerateNode(ModelPart& rTheVolumeModelPart, const Point& rCoordinates)`

# Calls

- [SetSolutionStepVariablesList](../../../../functions/includes/node/SetSolutionStepVariablesList.md)
- [pGetNodalSolutionStepVariablesList](../../../../functions/includes/model_part/pGetNodalSolutionStepVariablesList.md)