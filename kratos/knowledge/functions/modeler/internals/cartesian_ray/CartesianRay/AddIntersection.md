---
type: C++ Method
title: AddIntersection
resource: modeler/internals/cartesian_ray.h#L85-L99
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/modeler/internals/cartesian_ray/CartesianRay/ComputeTriangleRayIntersection
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/AddGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_ray/CartesianRay/CheckPassingThroughByExtraRays
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/modeler/test_voxel_mesh_modeler/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AddIntersection(TGeometryType const& rGeometry, double Tolerance)`

# Calls

- [ComputeTriangleRayIntersection](../../../../../functions/modeler/internals/cartesian_ray/CartesianRay/ComputeTriangleRayIntersection.md)

# Called by

- [AddGeometry](../../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/AddGeometry.md)
- [CheckPassingThroughByExtraRays](../../../../../functions/modeler/internals/cartesian_ray/CartesianRay/CheckPassingThroughByExtraRays.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/tests/cpp_tests/modeler/test_voxel_mesh_modeler/KRATOS_TEST_CASE_IN_SUITE.md)