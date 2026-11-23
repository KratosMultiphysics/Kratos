---
type: C++ Method
title: TransferObjectsToSonsNormalized
resource: spatial_containers/octree_binary_cell.h#L364-L394
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/spatial_containers/octree_binary_cell/OctreeBinaryCell/HasChildren
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/octree_binary_cell/OctreeBinaryCell/GetMinPointNormalized
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/octree_binary_cell/OctreeBinaryCell/GetMaxPointNormalized
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/spatial_containers/octree_binary/OctreeBinary/Constrain2To1
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void TransferObjectsToSonsNormalized()`

# Calls

- [HasChildren](../../../../functions/spatial_containers/octree_binary_cell/OctreeBinaryCell/HasChildren.md)
- [GetMinPointNormalized](../../../../functions/spatial_containers/octree_binary_cell/OctreeBinaryCell/GetMinPointNormalized.md)
- [GetMaxPointNormalized](../../../../functions/spatial_containers/octree_binary_cell/OctreeBinaryCell/GetMaxPointNormalized.md)

# Called by

- [Constrain2To1](../../../../functions/spatial_containers/octree_binary/OctreeBinary/Constrain2To1.md)