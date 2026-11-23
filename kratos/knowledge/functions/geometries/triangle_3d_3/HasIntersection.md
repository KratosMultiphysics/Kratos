---
type: C++ Function
title: HasIntersection
resource: geometries/triangle_3d_3.h#L669-L687
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/triangle_3d_3/LineTriangleOverlap
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/triangle_3d_3/TriangleTriangleOverlap
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/GeometryUtils/TriangleBoxOverlap
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool HasIntersection(const GeometryType& rThisGeometry) const override`

# Calls

- [LineTriangleOverlap](../../../functions/geometries/triangle_3d_3/LineTriangleOverlap.md)
- [TriangleTriangleOverlap](../../../functions/geometries/triangle_3d_3/TriangleTriangleOverlap.md)
- [TriangleBoxOverlap](../../../functions/utilities/geometry_utilities/GeometryUtils/TriangleBoxOverlap.md)