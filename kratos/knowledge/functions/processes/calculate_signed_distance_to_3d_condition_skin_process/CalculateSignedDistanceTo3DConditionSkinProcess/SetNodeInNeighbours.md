---
type: C++ Method
title: SetNodeInNeighbours
resource: processes/calculate_signed_distance_to_3d_condition_skin_process.h#L1350-L1372
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/spatial_containers/octree_binary_cell/OctreeBinaryCell/GetKey
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/octree_binary_cell/OctreeBinaryCell/GetLocalPosition
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SetNodeInNeighbours(CellType* pCell, int Position, CellNodeDataType* pNode)`

# Calls

- [GetKey](../../../../functions/spatial_containers/octree_binary_cell/OctreeBinaryCell/GetKey.md)
- [GetLocalPosition](../../../../functions/spatial_containers/octree_binary_cell/OctreeBinaryCell/GetLocalPosition.md)