---
type: C++ Function
title: SignedDistanceUnder
resource: modeler/internals/skin_intersection.h#L131-L154
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/internals/skin_intersection/TrimHigher
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static std::size_t SignedDistanceUnder( Vector& rNodalDistances, const GeometryType& rGeometry, const Point& rReferencePoint, const std::size_t Component)`

# Calls

- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)

# Called by

- [TrimHigher](../../../../functions/modeler/internals/skin_intersection/TrimHigher.md)