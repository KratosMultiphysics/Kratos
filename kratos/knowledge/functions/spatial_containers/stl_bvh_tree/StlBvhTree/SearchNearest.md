---
type: C++ Method
title: SearchNearest
resource: spatial_containers/stl_bvh_tree.cpp#L31-L164
visibility: private
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
---

# Signature

`StlBvhTree::NearestResult StlBvhTree::SearchNearest(const PointType& rCenter) const`

# Calls

- [AABBSphereDistance2](../../../../functions/spatial_containers/stl_bvh_tree/AABBSphereDistance2.md)
- [PackIndex](../../../../functions/spatial_containers/stl_bvh_tree/PackIndex.md)
- [PointTriangleDistance2](../../../../functions/spatial_containers/stl_bvh_tree/PointTriangleDistance2.md)