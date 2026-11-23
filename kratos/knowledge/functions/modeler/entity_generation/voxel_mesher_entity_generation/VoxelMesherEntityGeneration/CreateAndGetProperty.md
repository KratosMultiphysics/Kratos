---
type: C++ Method
title: CreateAndGetProperty
resource: modeler/entity_generation/voxel_mesher_entity_generation.cpp#L76-L82
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part/ModelPart/RecursivelyHasProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/CreateNewProperties
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

`Properties::Pointer VoxelMesherEntityGeneration::CreateAndGetProperty(ModelPart& rModelPart, std::size_t PropertyId) const`

# Calls

- [RecursivelyHasProperties](../../../../../functions/sources/model_part/ModelPart/RecursivelyHasProperties.md)
- [CreateNewProperties](../../../../../functions/sources/model_part/ModelPart/CreateNewProperties.md)

# Called by

- [Generate](../../../../../functions/modeler/entity_generation/generate_conditions_with_face_color/GenerateConditionsWithFaceColor/Generate.md)
- [Generate](../../../../../functions/modeler/entity_generation/generate_elements_with_cell_color/GenerateElementsWithCellColor/Generate.md)
- [Generate](../../../../../functions/modeler/entity_generation/generate_quadratic_elements_with_cell_color/GenerateQuadraticElementsWithCellColor/Generate.md)
- [Generate](../../../../../functions/modeler/entity_generation/generate_stl_intersection_with_cells/GenerateStlIntersectionWithCells/Generate.md)
- [Generate](../../../../../functions/modeler/entity_generation/generate_tetrahedral_elements_with_cell_color/GenerateTetrahedralElementsWithCellColor/Generate.md)
- [Generate](../../../../../functions/modeler/entity_generation/generate_triangular_conditions_with_face_color/GenerateTriangularConditionsWithFaceColor/Generate.md)