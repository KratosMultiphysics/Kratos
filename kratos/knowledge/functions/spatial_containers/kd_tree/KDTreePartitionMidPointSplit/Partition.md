---
type: C++ Method
title: Partition
resource: spatial_containers/kd_tree.h#L852-L878
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/spatial_containers/kd_tree/KDTreePartitionMidPointSplit/Split
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static IteratorType Partition( IteratorType PointsBegin, IteratorType PointsEnd, PointType const& HighPoint, PointType const& LowPoint, IndexType& rCuttingDimension, CoordinateType& rCuttingValue )`

# Calls

- [Split](../../../../functions/spatial_containers/kd_tree/KDTreePartitionMidPointSplit/Split.md)