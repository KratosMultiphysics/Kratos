---
type: C++ Method
title: Construct
resource: spatial_containers/kd_tree.h#L538-L576
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/spatial_containers/search_structure/PointerDistance
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static TreeNodeType* Construct( IteratorType PointsBegin, IteratorType PointsEnd, const PointType& HighPoint, const PointType& LowPoint, SizeType BucketSize )`

# Calls

- [PointerDistance](../../../../functions/spatial_containers/search_structure/PointerDistance.md)