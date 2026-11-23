---
type: C++ Method
title: Edge_Against_Tri_Edges
resource: geometries/triangle_2d_3.h#L1910-L1933
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/triangle_2d_3/Triangle2D3/Edge_Edge_Test
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/geometries/triangle_2d_3/Triangle2D3/coplanar_tri_tri
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool Edge_Against_Tri_Edges(const short& i0, const short& i1, const Point& V0, const Point& V1, const Point&U0, const Point&U1, const Point&U2) const`

# Calls

- [Edge_Edge_Test](../../../../functions/geometries/triangle_2d_3/Triangle2D3/Edge_Edge_Test.md)

# Called by

- [coplanar_tri_tri](../../../../functions/geometries/triangle_2d_3/Triangle2D3/coplanar_tri_tri.md)