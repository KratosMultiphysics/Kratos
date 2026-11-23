---
type: C++ Method
title: SetNumberOfDivisions
resource: modeler/voxel_mesh_generator_modeler.h#L137-L148
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/modeler/voxel_mesh_generator_modeler/CartesianNodalData/CartesianNodalData
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/entity_generation/voxel_mesher_entity_generation/VoxelMesherEntityGeneration/InitializeQuadraticData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/PreparingTheInternalDataStructure
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SetNumberOfDivisions(std::size_t XDivisions, std::size_t YDivisions, std::size_t ZDivisions)`

# Calls

- [CartesianNodalData](../../../../functions/modeler/voxel_mesh_generator_modeler/CartesianNodalData/CartesianNodalData.md)

# Called by

- [InitializeQuadraticData](../../../../functions/modeler/entity_generation/voxel_mesher_entity_generation/VoxelMesherEntityGeneration/InitializeQuadraticData.md)
- [PreparingTheInternalDataStructure](../../../../functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/PreparingTheInternalDataStructure.md)