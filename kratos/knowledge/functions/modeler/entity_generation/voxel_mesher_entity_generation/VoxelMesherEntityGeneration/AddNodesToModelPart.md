---
type: C++ Method
title: AddNodesToModelPart
resource: modeler/entity_generation/voxel_mesher_entity_generation.cpp#L111-L115
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/model_part_utils/AddNodesFromOrderedContainer
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/entity_generation/generate_conditions_with_face_color/GenerateConditionsWithFaceColor/Generate
    resolved_by: tree-sitter
    confidence: exact
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
  - target: functions/modeler/entity_generation/generate_triangular_conditions_with_face_color/GenerateTriangularConditionsWithFaceColor/Generate
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void VoxelMesherEntityGeneration::AddNodesToModelPart(ModelPart& rModelPart, ModelPart::NodesContainerType& rNewNodes) const`

# Calls

- [AddNodesFromOrderedContainer](../../../../../functions/utilities/model_part_utils/AddNodesFromOrderedContainer.md)

# Called by

- [Generate](../../../../../functions/modeler/entity_generation/generate_conditions_with_face_color/GenerateConditionsWithFaceColor/Generate.md)
- [Generate](../../../../../functions/modeler/entity_generation/generate_elements_with_cell_color/GenerateElementsWithCellColor/Generate.md)
- [Generate](../../../../../functions/modeler/entity_generation/generate_quadratic_elements_with_cell_color/GenerateQuadraticElementsWithCellColor/Generate.md)
- [Generate](../../../../../functions/modeler/entity_generation/generate_stl_intersection_with_cells/GenerateStlIntersectionWithCells/Generate.md)
- [Generate](../../../../../functions/modeler/entity_generation/generate_tetrahedral_elements_with_cell_color/GenerateTetrahedralElementsWithCellColor/Generate.md)
- [Generate](../../../../../functions/modeler/entity_generation/generate_triangular_conditions_with_face_color/GenerateTriangularConditionsWithFaceColor/Generate.md)