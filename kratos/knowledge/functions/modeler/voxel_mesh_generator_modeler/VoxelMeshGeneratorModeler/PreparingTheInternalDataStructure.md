---
type: C++ Method
title: PreparingTheInternalDataStructure
resource: modeler/voxel_mesh_generator_modeler.cpp#L141-L152
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/SetCoordinates
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/voxel_mesh_generator_modeler/CartesianData/SetNumberOfDivisions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_mesh_colors/ExtendBoundingBox
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/SetupModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void VoxelMeshGeneratorModeler::PreparingTheInternalDataStructure()`

# Calls

- [SetCoordinates](../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/SetCoordinates.md)
- [SetNumberOfDivisions](../../../../functions/modeler/voxel_mesh_generator_modeler/CartesianData/SetNumberOfDivisions.md)
- [ExtendBoundingBox](../../../../functions/modeler/internals/cartesian_mesh_colors/ExtendBoundingBox.md)

# Called by

- [SetupModelPart](../../../../functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/SetupModelPart.md)