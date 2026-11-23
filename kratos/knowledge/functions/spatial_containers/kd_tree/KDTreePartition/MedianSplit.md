---
type: C++ Method
title: MedianSplit
resource: spatial_containers/kd_tree.h#L480-L533
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/spatial_containers/kd_tree/KDTreePartition/Partition
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void MedianSplit( IteratorType PointsBegin, IteratorType PartitionPosition, IteratorType PointsEnd, IndexType CuttingDimension, CoordinateType& rCuttingValue )`

# Called by

- [Partition](../../../../functions/spatial_containers/kd_tree/KDTreePartition/Partition.md)