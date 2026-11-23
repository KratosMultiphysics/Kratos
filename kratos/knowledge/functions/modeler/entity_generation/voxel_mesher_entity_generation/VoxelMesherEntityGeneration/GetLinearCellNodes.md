---
type: C++ Method
title: GetLinearCellNodes
resource: modeler/entity_generation/voxel_mesher_entity_generation.cpp#L173-L189
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/modeler/entity_generation/generate_elements_with_cell_color/GenerateElementsWithCellColor/Generate
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/entity_generation/generate_tetrahedral_elements_with_cell_color/GenerateTetrahedralElementsWithCellColor/Generate
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/entity_generation/voxel_mesher_entity_generation/VoxelMesherEntityGeneration/GetQuadraticCellNodes
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void VoxelMesherEntityGeneration::GetLinearCellNodes( Element::NodesArrayType& rCellNodes, ModelPart& rTheVolumeModelPart, ModelPart::NodesContainerType& rThisNodes, const std::size_t I, const std::size_t J, const std::size_t K)`

# Called by

- [Generate](../../../../../functions/modeler/entity_generation/generate_elements_with_cell_color/GenerateElementsWithCellColor/Generate.md)
- [Generate](../../../../../functions/modeler/entity_generation/generate_tetrahedral_elements_with_cell_color/GenerateTetrahedralElementsWithCellColor/Generate.md)
- [GetQuadraticCellNodes](../../../../../functions/modeler/entity_generation/voxel_mesher_entity_generation/VoxelMesherEntityGeneration/GetQuadraticCellNodes.md)