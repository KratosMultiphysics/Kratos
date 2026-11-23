---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/modeler/test_voxel_mesh_modeler.cpp#L112-L152
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/cpp_tests/modeler/test_voxel_mesh_modeler/WriteCubeSkinMeshMdpaFileForVoxelModelerTest
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_ray/CartesianRay/AddIntersection
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_ray/CartesianRay/CollapseIntersectionPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_ray/CartesianRay/MarkIntersectedIntervals
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/SetCoordinates
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/InitializeRays
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/tetrahedra_mesh_edge_swapping_process/Edge/GetPoint1
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/tetrahedra_mesh_edge_swapping_process/Edge/GetPoint2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_mesh_colors/ExtendBoundingBox
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/SetAllColors
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/GetElementalColor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/CalculateElementalFaceColors
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(VoxelMeshGeneratorModelerKeyPlaneBySize, KratosCoreFastSuite)`

# Calls

- [WriteCubeSkinMeshMdpaFileForVoxelModelerTest](../../../../../functions/tests/cpp_tests/modeler/test_voxel_mesh_modeler/WriteCubeSkinMeshMdpaFileForVoxelModelerTest.md)
- [AddIntersection](../../../../../functions/modeler/internals/cartesian_ray/CartesianRay/AddIntersection.md)
- [CollapseIntersectionPoints](../../../../../functions/modeler/internals/cartesian_ray/CartesianRay/CollapseIntersectionPoints.md)
- [MarkIntersectedIntervals](../../../../../functions/modeler/internals/cartesian_ray/CartesianRay/MarkIntersectedIntervals.md)
- [SetCoordinates](../../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/SetCoordinates.md)
- [InitializeRays](../../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/InitializeRays.md)
- [GetPoint1](../../../../../functions/processes/tetrahedra_mesh_edge_swapping_process/Edge/GetPoint1.md)
- [GetPoint2](../../../../../functions/processes/tetrahedra_mesh_edge_swapping_process/Edge/GetPoint2.md)
- [ExtendBoundingBox](../../../../../functions/modeler/internals/cartesian_mesh_colors/ExtendBoundingBox.md)
- [GetGeometry](../../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [SetAllColors](../../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/SetAllColors.md)
- [GetElementalColor](../../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/GetElementalColor.md)
- [CalculateElementalFaceColors](../../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/CalculateElementalFaceColors.md)