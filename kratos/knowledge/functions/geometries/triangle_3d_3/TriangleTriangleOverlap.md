---
type: C++ Function
title: TriangleTriangleOverlap
resource: geometries/triangle_3d_3.h#L1974-L2056
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/plane_3d/Plane3D/CalculateSignedDistance
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/triangle_3d_3/AllSameSide
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/triangle_3d_3/GetMajorAxis
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/triangle_3d_3/ComputeIntervals
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/triangle_3d_3/CoplanarIntersectionCheck
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/geometries/triangle_3d_3/HasIntersection
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool TriangleTriangleOverlap( const Point& rPoint1, const Point& rPoint2, const Point& rPoint3) const`

# Calls

- [CalculateSignedDistance](../../../functions/geometries/plane_3d/Plane3D/CalculateSignedDistance.md)
- [AllSameSide](../../../functions/geometries/triangle_3d_3/AllSameSide.md)
- [GetMajorAxis](../../../functions/geometries/triangle_3d_3/GetMajorAxis.md)
- [ComputeIntervals](../../../functions/geometries/triangle_3d_3/ComputeIntervals.md)
- [CoplanarIntersectionCheck](../../../functions/geometries/triangle_3d_3/CoplanarIntersectionCheck.md)

# Called by

- [HasIntersection](../../../functions/geometries/triangle_3d_3/HasIntersection.md)