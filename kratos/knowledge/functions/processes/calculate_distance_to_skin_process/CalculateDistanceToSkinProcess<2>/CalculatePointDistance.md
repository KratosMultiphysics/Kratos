---
type: C++ Method
title: CalculatePointDistance
resource: processes/calculate_distance_to_skin_process.cpp#L277-L285
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/geometry_utilities/GeometryUtils/PointDistanceToLineSegment3D
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`double inline CalculateDistanceToSkinProcess<2>::CalculatePointDistance( const Element::GeometryType &rIntObjGeom, const Point &rDistancePoint)`

# Calls

- [PointDistanceToLineSegment3D](../../../../functions/utilities/geometry_utilities/GeometryUtils/PointDistanceToLineSegment3D.md)