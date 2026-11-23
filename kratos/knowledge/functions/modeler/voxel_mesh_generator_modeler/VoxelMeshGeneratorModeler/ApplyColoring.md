---
type: C++ Method
title: ApplyColoring
resource: modeler/voxel_mesh_generator_modeler.cpp#L154-L177
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/timer/Timer/Start
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/SetAllColors
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/timer/Timer/Stop
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/SetupModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void VoxelMeshGeneratorModeler::ApplyColoring(Parameters ColoringParameters, int OutsideColor)`

# Calls

- [Start](../../../../functions/utilities/timer/Timer/Start.md)
- [SetAllColors](../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/SetAllColors.md)
- [AddInt](../../../../functions/sources/kratos_parameters/Parameters/AddInt.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [Stop](../../../../functions/utilities/timer/Timer/Stop.md)

# Called by

- [SetupModelPart](../../../../functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/SetupModelPart.md)