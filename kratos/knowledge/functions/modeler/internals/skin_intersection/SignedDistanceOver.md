---
type: C++ Function
title: SignedDistanceOver
resource: modeler/internals/skin_intersection.h#L159-L182
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/internals/skin_intersection/TrimLower
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static std::size_t SignedDistanceOver( Vector& rNodalDistances, const GeometryType& rGeometry, const Point& rReferencePoint, const std::size_t Component)`

# Calls

- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)

# Called by

- [TrimLower](../../../../functions/modeler/internals/skin_intersection/TrimLower.md)