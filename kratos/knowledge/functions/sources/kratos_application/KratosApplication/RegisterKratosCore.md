---
type: C++ Method
title: RegisterKratosCore
resource: sources/kratos_application.cpp#L153-L361
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/variables/KratosApplication/RegisterVariables
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/factories/standard_preconditioner_factory/RegisterPreconditioners
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/coloring/voxel_mesher_coloring_factory/RegisterVoxelMesherColoring
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/key_plane_generation/key_plane_generation_factory/RegisterVoxelMesherKeyPlaneGeneration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/entity_generation/voxel_mesher_entity_generation_factory/RegisterVoxelMesherEntityGeneration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/operation/voxel_mesher_operation_factory/RegisterVoxelMesherOperation
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/includes/kratos_application/Register
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void KratosApplication::RegisterKratosCore()`

# Calls

- [RegisterVariables](../../../../functions/sources/variables/KratosApplication/RegisterVariables.md)
- [RegisterPreconditioners](../../../../functions/factories/standard_preconditioner_factory/RegisterPreconditioners.md)
- [RegisterVoxelMesherColoring](../../../../functions/modeler/coloring/voxel_mesher_coloring_factory/RegisterVoxelMesherColoring.md)
- [RegisterVoxelMesherKeyPlaneGeneration](../../../../functions/modeler/key_plane_generation/key_plane_generation_factory/RegisterVoxelMesherKeyPlaneGeneration.md)
- [RegisterVoxelMesherEntityGeneration](../../../../functions/modeler/entity_generation/voxel_mesher_entity_generation_factory/RegisterVoxelMesherEntityGeneration.md)
- [RegisterVoxelMesherOperation](../../../../functions/modeler/operation/voxel_mesher_operation_factory/RegisterVoxelMesherOperation.md)

# Called by

- [Register](../../../../functions/includes/kratos_application/Register.md)