---
type: C++ Method
title: OCTreePartition
resource: spatial_containers/octree.h#L130-L216
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/spatial_containers/search_structure/PointerDistance
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/kd_tree/KDTreePartitionAverageSplit/AverageSplit
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`OCTreePartition(IteratorType PointsBegin, IteratorType PointsEnd, PointType const& MinPoint, PointType const& MaxPoint, SizeType BucketSize = 1)`

# Calls

- [PointerDistance](../../../../functions/spatial_containers/search_structure/PointerDistance.md)
- [AverageSplit](../../../../functions/spatial_containers/kd_tree/KDTreePartitionAverageSplit/AverageSplit.md)