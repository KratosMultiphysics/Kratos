---
type: C++ Function
title: ApplyInRadius
resource: spatial_containers/stl_bvh_tree.h#L103-L154
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/spatial_containers/stl_bvh_tree/AABBSphereDistance2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/stl_bvh_tree/PackIndex
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/stl_bvh_tree/PointTriangleDistance2
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/spatial_containers/test_stl_bvh_tree/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ApplyInRadius( const PointType& rCenter, double Radius, TFunction&& rFunction) const`

# Calls

- [AABBSphereDistance2](../../../functions/spatial_containers/stl_bvh_tree/AABBSphereDistance2.md)
- [PackIndex](../../../functions/spatial_containers/stl_bvh_tree/PackIndex.md)
- [PointTriangleDistance2](../../../functions/spatial_containers/stl_bvh_tree/PointTriangleDistance2.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/spatial_containers/test_stl_bvh_tree/KRATOS_TEST_CASE_IN_SUITE.md)