---
type: C++ Function
title: TrimLower
resource: modeler/internals/skin_intersection.h#L214-L240
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/modeler/internals/skin_intersection/SignedDistanceOver
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetPositiveSubdivisions
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static std::vector<GeometryPointerType> TrimLower( const std::vector<GeometryPointerType>& rGeometries, const Point& rHighPoint, std::size_t Component)`

# Calls

- [SignedDistanceOver](../../../../functions/modeler/internals/skin_intersection/SignedDistanceOver.md)
- [GetPositiveSubdivisions](../../../../functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetPositiveSubdivisions.md)