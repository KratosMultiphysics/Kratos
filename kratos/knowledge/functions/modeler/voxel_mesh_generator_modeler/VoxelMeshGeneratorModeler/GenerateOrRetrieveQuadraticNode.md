---
type: C++ Method
title: GenerateOrRetrieveQuadraticNode
resource: modeler/voxel_mesh_generator_modeler.cpp#L241-L260
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/modeler/voxel_mesh_generator_modeler/CartesianNodalData/IsCreated
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/voxel_mesh_generator_modeler/CartesianNodalData/pSetNode
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`Node::Pointer VoxelMeshGeneratorModeler::GenerateOrRetrieveQuadraticNode( ModelPart& rTheVolumeModelPart, ModelPart::NodesContainerType& rThisNodes, const std::size_t I, const std::size_t J, const std::size_t K )`

# Calls

- [IsCreated](../../../../functions/modeler/voxel_mesh_generator_modeler/CartesianNodalData/IsCreated.md)
- [pSetNode](../../../../functions/modeler/voxel_mesh_generator_modeler/CartesianNodalData/pSetNode.md)