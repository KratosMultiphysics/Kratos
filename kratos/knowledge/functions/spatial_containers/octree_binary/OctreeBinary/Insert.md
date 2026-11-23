---
type: C++ Method
title: Insert
resource: spatial_containers/octree_binary.h#L281-L297
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/spatial_containers/octree_binary/OctreeBinary/CalcKeyNormalized
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/pointer_vector_set/PointerVectorSet/pop_back
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/octree_binary_cell/OctreeBinaryCell/HasChildren
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/octree_binary_cell/OctreeBinaryCell/GetMinPointNormalized
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/octree_binary_cell/OctreeBinaryCell/GetMaxPointNormalized
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Insert(coordinate_type* point)`

# Calls

- [CalcKeyNormalized](../../../../functions/spatial_containers/octree_binary/OctreeBinary/CalcKeyNormalized.md)
- [pop_back](../../../../functions/containers/pointer_vector_set/PointerVectorSet/pop_back.md)
- [HasChildren](../../../../functions/spatial_containers/octree_binary_cell/OctreeBinaryCell/HasChildren.md)
- [GetMinPointNormalized](../../../../functions/spatial_containers/octree_binary_cell/OctreeBinaryCell/GetMinPointNormalized.md)
- [GetMaxPointNormalized](../../../../functions/spatial_containers/octree_binary_cell/OctreeBinaryCell/GetMaxPointNormalized.md)