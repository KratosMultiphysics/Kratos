---
type: C++ Method
title: Generate
resource: modeler/entity_generation/generate_tetrahedral_elements_with_cell_color.cpp#L34-L66
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/modeler/entity_generation/voxel_mesher_entity_generation/VoxelMesherEntityGeneration/GetStartElementId
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/entity_generation/voxel_mesher_entity_generation/VoxelMesherEntityGeneration/CreateAndGetProperty
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/entity_generation/voxel_mesher_entity_generation/VoxelMesherEntityGeneration/GetNumberOfCells
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/GetElementalColor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/entity_generation/voxel_mesher_entity_generation/VoxelMesherEntityGeneration/GetLinearCellNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/entity_generation/generate_tetrahedral_elements_with_cell_color/GenerateTetrahedralElementsWithCellColor/CreateTetrahedraInCell
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/entity_generation/voxel_mesher_entity_generation/VoxelMesherEntityGeneration/AddNodesToModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void GenerateTetrahedralElementsWithCellColor::Generate(ModelPart& rModelPart, Parameters parameters)`

# Calls

- [GetStartElementId](../../../../../functions/modeler/entity_generation/voxel_mesher_entity_generation/VoxelMesherEntityGeneration/GetStartElementId.md)
- [GetInt](../../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [CreateAndGetProperty](../../../../../functions/modeler/entity_generation/voxel_mesher_entity_generation/VoxelMesherEntityGeneration/CreateAndGetProperty.md)
- [GetString](../../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [GetNumberOfCells](../../../../../functions/modeler/entity_generation/voxel_mesher_entity_generation/VoxelMesherEntityGeneration/GetNumberOfCells.md)
- [GetElementalColor](../../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/GetElementalColor.md)
- [GetLinearCellNodes](../../../../../functions/modeler/entity_generation/voxel_mesher_entity_generation/VoxelMesherEntityGeneration/GetLinearCellNodes.md)
- [CreateTetrahedraInCell](../../../../../functions/modeler/entity_generation/generate_tetrahedral_elements_with_cell_color/GenerateTetrahedralElementsWithCellColor/CreateTetrahedraInCell.md)
- [AddNodesToModelPart](../../../../../functions/modeler/entity_generation/voxel_mesher_entity_generation/VoxelMesherEntityGeneration/AddNodesToModelPart.md)