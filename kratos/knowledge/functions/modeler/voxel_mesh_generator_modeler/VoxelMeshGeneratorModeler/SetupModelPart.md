---
type: C++ Method
title: SetupModelPart
resource: modeler/voxel_mesh_generator_modeler.cpp#L82-L117
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/PreparingTheInternalDataStructure
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/ApplyColoring
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/GenerateEntities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/ApplyOperations
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/WriteParaViewVTR
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/PrintAllVoxelColorsToAsciiFile
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void VoxelMeshGeneratorModeler::SetupModelPart()`

# Calls

- [PreparingTheInternalDataStructure](../../../../functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/PreparingTheInternalDataStructure.md)
- [ApplyColoring](../../../../functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/ApplyColoring.md)
- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [GenerateEntities](../../../../functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/GenerateEntities.md)
- [ApplyOperations](../../../../functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/ApplyOperations.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [WriteParaViewVTR](../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/WriteParaViewVTR.md)
- [PrintAllVoxelColorsToAsciiFile](../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/PrintAllVoxelColorsToAsciiFile.md)