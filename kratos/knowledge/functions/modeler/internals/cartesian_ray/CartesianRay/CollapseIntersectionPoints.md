---
type: C++ Method
title: CollapseIntersectionPoints
resource: modeler/internals/cartesian_ray.h#L101-L136
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/modeler/internals/cartesian_ray/CartesianRay/CheckPassingThroughByExtraRays
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/CalculateNodalRayColors
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/CalculateElementalRayColors
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/CalculateElementalFaceColors
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/modeler/test_voxel_mesh_modeler/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CollapseIntersectionPoints(double Tolerance)`

# Calls

- [CheckPassingThroughByExtraRays](../../../../../functions/modeler/internals/cartesian_ray/CartesianRay/CheckPassingThroughByExtraRays.md)

# Called by

- [CalculateNodalRayColors](../../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/CalculateNodalRayColors.md)
- [CalculateElementalRayColors](../../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/CalculateElementalRayColors.md)
- [CalculateElementalFaceColors](../../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/CalculateElementalFaceColors.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/tests/cpp_tests/modeler/test_voxel_mesh_modeler/KRATOS_TEST_CASE_IN_SUITE.md)