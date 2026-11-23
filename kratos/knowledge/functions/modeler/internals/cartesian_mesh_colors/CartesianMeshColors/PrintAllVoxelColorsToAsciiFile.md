---
type: C++ Method
title: PrintAllVoxelColorsToAsciiFile
resource: modeler/internals/cartesian_mesh_colors.h#L774-L799
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/GetElementalColor
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/SetupModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void PrintAllVoxelColorsToAsciiFile(const std::string& filename)`

# Calls

- [GetElementalColor](../../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/GetElementalColor.md)

# Called by

- [SetupModelPart](../../../../../functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/SetupModelPart.md)