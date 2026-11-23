---
type: C++ Method
title: BuildNodeInternal
resource: spatial_containers/stl_bvh_tree.cpp#L170-L282
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/spatial_containers/stl_bvh_tree/PackRef
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/pointer_vector/emplace_back
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/spatial_containers/stl_bvh_tree/StlBvhTree/Build
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`StlBvhTree::BuildResult StlBvhTree::BuildNodeInternal(int PackBegin, int PackEnd)`

# Calls

- [PackRef](../../../../functions/spatial_containers/stl_bvh_tree/PackRef.md)
- [emplace_back](../../../../functions/containers/pointer_vector/emplace_back.md)

# Called by

- [Build](../../../../functions/spatial_containers/stl_bvh_tree/StlBvhTree/Build.md)