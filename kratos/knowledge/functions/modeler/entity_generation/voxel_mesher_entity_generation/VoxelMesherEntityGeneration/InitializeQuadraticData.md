---
type: C++ Method
title: InitializeQuadraticData
resource: modeler/entity_generation/voxel_mesher_entity_generation.cpp#L239-L244
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/modeler/voxel_mesh_generator_modeler/CartesianData/SetNumberOfDivisions
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/entity_generation/generate_quadratic_elements_with_cell_color/GenerateQuadraticElementsWithCellColor/Generate
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void VoxelMesherEntityGeneration::InitializeQuadraticData()`

# Calls

- [SetNumberOfDivisions](../../../../../functions/modeler/voxel_mesh_generator_modeler/CartesianData/SetNumberOfDivisions.md)

# Called by

- [Generate](../../../../../functions/modeler/entity_generation/generate_quadratic_elements_with_cell_color/GenerateQuadraticElementsWithCellColor/Generate.md)