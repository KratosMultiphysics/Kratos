---
type: C++ Method
title: SearchInRadius
resource: spatial_containers/octree.h#L258-L273
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/spatial_containers/octree/OCTreePartition/DistanceToPartitions
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SearchInRadius(PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results, DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults) override`

# Calls

- [DistanceToPartitions](../../../../functions/spatial_containers/octree/OCTreePartition/DistanceToPartitions.md)