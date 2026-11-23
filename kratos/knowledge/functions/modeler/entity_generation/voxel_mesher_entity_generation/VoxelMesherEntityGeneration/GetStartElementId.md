---
type: C++ Method
title: GetStartElementId
resource: modeler/entity_generation/voxel_mesher_entity_generation.cpp#L124-L127
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/modeler/entity_generation/generate_elements_with_cell_color/GenerateElementsWithCellColor/Generate
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/entity_generation/generate_quadratic_elements_with_cell_color/GenerateQuadraticElementsWithCellColor/Generate
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/entity_generation/generate_stl_intersection_with_cells/GenerateStlIntersectionWithCells/Generate
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/entity_generation/generate_tetrahedral_elements_with_cell_color/GenerateTetrahedralElementsWithCellColor/Generate
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::size_t VoxelMesherEntityGeneration::GetStartElementId() const`

# Called by

- [Generate](../../../../../functions/modeler/entity_generation/generate_elements_with_cell_color/GenerateElementsWithCellColor/Generate.md)
- [Generate](../../../../../functions/modeler/entity_generation/generate_quadratic_elements_with_cell_color/GenerateQuadraticElementsWithCellColor/Generate.md)
- [Generate](../../../../../functions/modeler/entity_generation/generate_stl_intersection_with_cells/GenerateStlIntersectionWithCells/Generate.md)
- [Generate](../../../../../functions/modeler/entity_generation/generate_tetrahedral_elements_with_cell_color/GenerateTetrahedralElementsWithCellColor/Generate.md)