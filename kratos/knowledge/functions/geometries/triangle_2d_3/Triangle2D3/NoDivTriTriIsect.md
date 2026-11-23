---
type: C++ Method
title: NoDivTriTriIsect
resource: geometries/triangle_2d_3.h#L1625-L1753
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/triangle_2d_3/Triangle2D3/New_Compute_Intervals
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/triangle_2d_3/Triangle2D3/coplanar_tri_tri
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/geometries/triangle_2d_3/Triangle2D3/HasIntersection
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool NoDivTriTriIsect( const Point& V0, const Point& V1, const Point& V2, const Point& U0, const Point& U1, const Point& U2) const`

# Calls

- [New_Compute_Intervals](../../../../functions/geometries/triangle_2d_3/Triangle2D3/New_Compute_Intervals.md)
- [coplanar_tri_tri](../../../../functions/geometries/triangle_2d_3/Triangle2D3/coplanar_tri_tri.md)

# Called by

- [HasIntersection](../../../../functions/geometries/triangle_2d_3/Triangle2D3/HasIntersection.md)