---
type: C++ Method
title: DistanceToPartitions
resource: spatial_containers/octree.h#L307-L330
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/spatial_containers/octree/OCTreePartition/SearchNearestPoint
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/octree/OCTreePartition/SearchInRadius
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void DistanceToPartitions(IndexType ContainingChildIndex, PointType const& rThisPoint, CoordinateType rDistances[]) const`

# Called by

- [SearchNearestPoint](../../../../functions/spatial_containers/octree/OCTreePartition/SearchNearestPoint.md)
- [SearchInRadius](../../../../functions/spatial_containers/octree/OCTreePartition/SearchInRadius.md)