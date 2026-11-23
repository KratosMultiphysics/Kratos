---
type: C++ Module
title: cartesian_mesh_colors
resource: modeler/internals/cartesian_mesh_colors.h#L1-L813
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/string
    resolved_by: tree-sitter
    confidence: exact
  - target: external/iostream
    resolved_by: tree-sitter
    confidence: exact
  - target: external/algorithm
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-model-part-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/modeler-internals-cartesian-ray-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-parallel-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [CartesianMeshColors](../../../classes/modeler/internals/cartesian_mesh_colors/CartesianMeshColors.md)
- [ElementIndex](../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/ElementIndex.md)
- [CartesianMeshColors](../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/CartesianMeshColors.md)
- [ExtendBoundingBox](../../../functions/modeler/internals/cartesian_mesh_colors/ExtendBoundingBox.md)
- [SetCoordinates](../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/SetCoordinates.md)
- [SetAllColors](../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/SetAllColors.md)
- [GetElementalColor](../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/GetElementalColor.md)
- [GetPoint](../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/GetPoint.md)
- [GetPoint](../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/GetPoint-2.md)
- [GetCenterOfElement](../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/GetCenterOfElement.md)
- [GetGridSize](../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/GetGridSize.md)
- [InitializeRays](../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/InitializeRays.md)
- [AddGeometry](../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/AddGeometry.md)
- [CalculateNodalRayColors](../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/CalculateNodalRayColors.md)
- [CalculateElementalRayColors](../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/CalculateElementalRayColors.md)
- [CalculateElementalFaceColors](../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/CalculateElementalFaceColors.md)
- [CalculateElementalFaceColorsBetweenColors](../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/CalculateElementalFaceColorsBetweenColors.md)
- [CalculateOuterMinMaxNodePositions](../../../functions/modeler/internals/cartesian_mesh_colors/CalculateOuterMinMaxNodePositions.md)
- [CalculateMinMaxNodePositions](../../../functions/modeler/internals/cartesian_mesh_colors/CalculateMinMaxNodePositions.md)
- [CalculateMinMaxPoints](../../../functions/modeler/internals/cartesian_mesh_colors/CalculateMinMaxPoints.md)
- [CalculateMinMaxCenterOfElementPositions](../../../functions/modeler/internals/cartesian_mesh_colors/CalculateMinMaxCenterOfElementPositions.md)
- [CalculateNodePosition](../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/CalculateNodePosition.md)
- [CalculateCenterOfElementPosition](../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/CalculateCenterOfElementPosition.md)
- [WriteParaViewVTR](../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/WriteParaViewVTR.md)
- [PrintAllVoxelColorsToAsciiFile](../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/PrintAllVoxelColorsToAsciiFile.md)

# Imports

- `<string>`
- `<iostream>`
- `<algorithm>`
- `"processes/process.h"`
- `"includes/model_part.h"`
- `"modeler/internals/cartesian_ray.h"`
- `"utilities/parallel_utilities.h"`

# Member of

- [KratosMultiphysics](../../../packages/KratosMultiphysics.md)