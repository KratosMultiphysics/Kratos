---
type: C++ Function
title: CoplanarIntersectionCheck
resource: geometries/triangle_3d_3.h#L2120-L2163
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/triangle_3d_3/EdgeToTriangleEdgesCheck
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/geometries/triangle_3d_3/TriangleTriangleOverlap
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool CoplanarIntersectionCheck( const array_1d<double,3>& N, const Point& rPoint1, const Point& rPoint2, const Point& rPoint3) const`

# Calls

- [EdgeToTriangleEdgesCheck](../../../functions/geometries/triangle_3d_3/EdgeToTriangleEdgesCheck.md)

# Called by

- [TriangleTriangleOverlap](../../../functions/geometries/triangle_3d_3/TriangleTriangleOverlap.md)