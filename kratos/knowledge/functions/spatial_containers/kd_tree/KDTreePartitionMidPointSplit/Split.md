---
type: C++ Method
title: Split
resource: spatial_containers/kd_tree.h#L923-L942
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/spatial_containers/kd_tree/KDTreePartitionMidPointSplit/Partition
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static IteratorType Split( IteratorType PointsBegin, IteratorType PointsEnd, IndexType& CuttingDimension, CoordinateType& rCuttingValue )`

# Called by

- [Partition](../../../../functions/spatial_containers/kd_tree/KDTreePartitionMidPointSplit/Partition.md)