---
type: C++ Method
title: SearchNearestInBoxInner
resource: spatial_containers/bins_static.h#L747-L752
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/spatial_containers/bins_static/Bins/SearchNearestInnerInRange
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SearchNearestInBoxInner( PointerType& ThisPoint, PointerType& ResultPoint, CoordinateType& ResultDistance, SearchStructure<IndexType,std::size_t,CoordinateType,IteratorType,IteratorIteratorType,1>& Box, bool& Found )`

# Calls

- [SearchNearestInnerInRange](../../../../functions/spatial_containers/bins_static/Bins/SearchNearestInnerInRange.md)