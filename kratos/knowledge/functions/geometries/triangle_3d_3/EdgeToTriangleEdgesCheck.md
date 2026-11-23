---
type: C++ Function
title: EdgeToTriangleEdgesCheck
resource: geometries/triangle_3d_3.h#L2165-L2188
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/triangle_3d_3/EdgeToEdgeIntersectionCheck
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/geometries/triangle_3d_3/CoplanarIntersectionCheck
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool EdgeToTriangleEdgesCheck( const int& i0, const int& i1, const Point& V0, const Point& V1, const Point&U0, const Point&U1, const Point&U2) const`

# Calls

- [EdgeToEdgeIntersectionCheck](../../../functions/geometries/triangle_3d_3/EdgeToEdgeIntersectionCheck.md)

# Called by

- [CoplanarIntersectionCheck](../../../functions/geometries/triangle_3d_3/CoplanarIntersectionCheck.md)