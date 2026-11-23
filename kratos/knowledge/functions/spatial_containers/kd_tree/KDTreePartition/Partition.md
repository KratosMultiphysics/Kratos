---
type: C++ Method
title: Partition
resource: spatial_containers/kd_tree.h#L419-L434
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/spatial_containers/search_structure/PointerDistance
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/kd_tree/KDTreePartition/MedianSplit
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static IteratorType Partition( IteratorType PointsBegin, IteratorType PointsEnd, IndexType& rCuttingDimension, CoordinateType& rCuttingValue )`

# Calls

- [PointerDistance](../../../../functions/spatial_containers/search_structure/PointerDistance.md)
- [MedianSplit](../../../../functions/spatial_containers/kd_tree/KDTreePartition/MedianSplit.md)