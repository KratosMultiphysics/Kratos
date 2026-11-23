---
type: C++ Method
title: GenerateEntities
resource: modeler/voxel_mesh_generator_modeler.cpp#L192-L212
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/timer/Timer/Start
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

`void VoxelMeshGeneratorModeler::GenerateEntities( ModelPart& rTheVolumeModelPart, Parameters EntityGeneratorParameters )`

# Calls

- [Start](../../../../functions/utilities/timer/Timer/Start.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [Stop](../../../../functions/utilities/timer/Timer/Stop.md)

# Called by

- [SetupModelPart](../../../../functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/SetupModelPart.md)