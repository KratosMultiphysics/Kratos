---
type: C++ Module
title: voxel_mesh_generator_modeler
resource: modeler/voxel_mesh_generator_modeler.cpp#L1-L321
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/includes-mesh-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-element-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-model-part-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-bounding-box-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-point-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-parallel-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-timer-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/spatial-containers-geometrical-objects-bins-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-model-part-utils-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/voxel-mesh-generator-modeler-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/coloring-voxel-mesher-coloring-factory-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/key-plane-generation-key-plane-generation-factory-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/entity-generation-voxel-mesher-entity-generation-factory-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/operation-voxel-mesher-operation-factory-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [VoxelMeshGeneratorModeler](../../functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/VoxelMeshGeneratorModeler.md)
- [VoxelMeshGeneratorModeler](../../functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/VoxelMeshGeneratorModeler-2.md)
- [GetDefaultParameters](../../functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/GetDefaultParameters.md)
- [GetDefaultKeyPlaneGeneratorParameters](../../functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/GetDefaultKeyPlaneGeneratorParameters.md)
- [SetupModelPart](../../functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/SetupModelPart.md)
- [PreparingTheInternalDataStructure](../../functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/PreparingTheInternalDataStructure.md)
- [ApplyColoring](../../functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/ApplyColoring.md)
- [SetStartIds](../../functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/SetStartIds.md)
- [GenerateEntities](../../functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/GenerateEntities.md)
- [GenerateOrRetrieveNode](../../functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/GenerateOrRetrieveNode.md)
- [GenerateOrRetrieveQuadraticNode](../../functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/GenerateOrRetrieveQuadraticNode.md)
- [GenerateNode](../../functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/GenerateNode.md)
- [ApplyOperations](../../functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/ApplyOperations.md)
- [IntersectionsGenerated](../../functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/IntersectionsGenerated.md)

# Imports

- `"includes/mesh.h"`
- `"includes/element.h"`
- `"includes/model_part.h"`
- `"geometries/bounding_box.h"`
- `"geometries/point.h"`
- `"utilities/parallel_utilities.h"`
- `"utilities/timer.h"`
- `"spatial_containers/geometrical_objects_bins.h"`
- `"utilities/model_part_utils.h"`
- `"voxel_mesh_generator_modeler.h"`
- `"coloring/voxel_mesher_coloring_factory.h"`
- `"key_plane_generation/key_plane_generation_factory.h"`
- `"entity_generation/voxel_mesher_entity_generation_factory.h"`
- `"operation/voxel_mesher_operation_factory.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)