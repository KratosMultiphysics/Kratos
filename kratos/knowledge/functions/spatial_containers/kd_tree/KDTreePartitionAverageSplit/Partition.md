---
type: C++ Method
title: Partition
resource: spatial_containers/kd_tree.h#L657-L672
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/spatial_containers/kd_tree/KDTreePartitionAverageSplit/AverageSplit
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static IteratorType Partition( IteratorType PointsBegin, IteratorType PointsEnd, IndexType& rCuttingDimension, CoordinateType& rCuttingValue )`

# Calls

- [AverageSplit](../../../../functions/spatial_containers/kd_tree/KDTreePartitionAverageSplit/AverageSplit.md)