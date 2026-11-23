---
type: C++ Method
title: AverageSplit
resource: spatial_containers/kd_tree.h#L725-L745
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/spatial_containers/kd_tree/KDTreePartitionAverageSplit/Partition
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/octree/OCTreePartition/OCTreePartition
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void AverageSplit( IteratorType PointsBegin, IteratorType& PartitionPosition, IteratorType PointsEnd, IndexType& CuttingDimension, CoordinateType& rCuttingValue )`

# Called by

- [Partition](../../../../functions/spatial_containers/kd_tree/KDTreePartitionAverageSplit/Partition.md)
- [OCTreePartition](../../../../functions/spatial_containers/octree/OCTreePartition/OCTreePartition.md)