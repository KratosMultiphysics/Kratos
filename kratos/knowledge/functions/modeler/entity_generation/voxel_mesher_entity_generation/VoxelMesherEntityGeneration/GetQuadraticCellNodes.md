---
type: C++ Method
title: GetQuadraticCellNodes
resource: modeler/entity_generation/voxel_mesher_entity_generation.cpp#L192-L236
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/modeler/entity_generation/voxel_mesher_entity_generation/VoxelMesherEntityGeneration/GetLinearCellNodes
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/entity_generation/generate_quadratic_elements_with_cell_color/GenerateQuadraticElementsWithCellColor/Generate
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void VoxelMesherEntityGeneration::GetQuadraticCellNodes( Element::NodesArrayType& rCellNodes, ModelPart& rTheVolumeModelPart, ModelPart::NodesContainerType& rThisNodes, const std::size_t I, const std::size_t J, const std::size_t K)`

# Calls

- [GetLinearCellNodes](../../../../../functions/modeler/entity_generation/voxel_mesher_entity_generation/VoxelMesherEntityGeneration/GetLinearCellNodes.md)

# Called by

- [Generate](../../../../../functions/modeler/entity_generation/generate_quadratic_elements_with_cell_color/GenerateQuadraticElementsWithCellColor/Generate.md)