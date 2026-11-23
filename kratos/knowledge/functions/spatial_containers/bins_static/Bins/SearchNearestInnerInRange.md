---
type: C++ Method
title: SearchNearestInnerInRange
resource: spatial_containers/bins_static.h#L776-L790
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/spatial_containers/bins_static/Bins/SearchNearestInBoxInner
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SearchNearestInnerInRange( const IteratorType& RangeBegin, const IteratorType& RangeEnd, PointerType& ThisPoint, PointerType& Result, CoordinateType& Distance, bool& Found )`

# Called by

- [SearchNearestInBoxInner](../../../../functions/spatial_containers/bins_static/Bins/SearchNearestInBoxInner.md)