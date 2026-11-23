---
type: C++ Function
title: AABBSphereDistance2
resource: spatial_containers/stl_bvh_tree.h#L358-L367
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

`static float AABBSphereDistance2( float cx, float cy, float cz, float lx, float ly, float lz, float ux, float uy, float uz)`

# Called by

- [SearchNearest](../../../functions/spatial_containers/stl_bvh_tree/StlBvhTree/SearchNearest.md)
- [ApplyInRadius](../../../functions/spatial_containers/stl_bvh_tree/ApplyInRadius.md)