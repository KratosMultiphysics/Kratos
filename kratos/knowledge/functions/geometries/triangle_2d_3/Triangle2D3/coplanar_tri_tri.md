---
type: C++ Method
title: coplanar_tri_tri
resource: geometries/triangle_2d_3.h#L1847-L1905
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/triangle_2d_3/Triangle2D3/Edge_Against_Tri_Edges
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/triangle_2d_3/Triangle2D3/Point_In_Tri
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/geometries/triangle_2d_3/Triangle2D3/NoDivTriTriIsect
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool coplanar_tri_tri( const array_1d<double, 3>& N, const Point& V0, const Point& V1, const Point& V2, const Point& U0, const Point& U1, const Point& U2) const`

# Calls

- [Edge_Against_Tri_Edges](../../../../functions/geometries/triangle_2d_3/Triangle2D3/Edge_Against_Tri_Edges.md)
- [Point_In_Tri](../../../../functions/geometries/triangle_2d_3/Triangle2D3/Point_In_Tri.md)

# Called by

- [NoDivTriTriIsect](../../../../functions/geometries/triangle_2d_3/Triangle2D3/NoDivTriTriIsect.md)