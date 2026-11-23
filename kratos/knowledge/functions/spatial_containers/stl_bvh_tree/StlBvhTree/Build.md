---
type: C++ Method
title: Build
resource: spatial_containers/stl_bvh_tree.cpp#L284-L381
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/spatial_containers/stl_bvh_tree/EncodeMorton3
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/pointer_vector/emplace_back
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/stl_bvh_tree/PackRef
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/stl_bvh_tree/StlBvhTree/BuildNodeInternal
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void StlBvhTree::Build(std::vector<TriangleData>& rRawTris)`

# Calls

- [EncodeMorton3](../../../../functions/spatial_containers/stl_bvh_tree/EncodeMorton3.md)
- [emplace_back](../../../../functions/containers/pointer_vector/emplace_back.md)
- [PackRef](../../../../functions/spatial_containers/stl_bvh_tree/PackRef.md)
- [BuildNodeInternal](../../../../functions/spatial_containers/stl_bvh_tree/StlBvhTree/BuildNodeInternal.md)