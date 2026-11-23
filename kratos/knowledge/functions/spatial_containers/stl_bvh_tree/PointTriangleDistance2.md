---
type: C++ Function
title: PointTriangleDistance2
resource: spatial_containers/stl_bvh_tree.h#L271-L355
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/spatial_containers/stl_bvh_tree/StlBvhTree/SearchNearest
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/stl_bvh_tree/ApplyInRadius
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static float PointTriangleDistance2( float px, float py, float pz, const float v0[3], const float e1[3], const float e2[3])`

# Called by

- [SearchNearest](../../../functions/spatial_containers/stl_bvh_tree/StlBvhTree/SearchNearest.md)
- [ApplyInRadius](../../../functions/spatial_containers/stl_bvh_tree/ApplyInRadius.md)