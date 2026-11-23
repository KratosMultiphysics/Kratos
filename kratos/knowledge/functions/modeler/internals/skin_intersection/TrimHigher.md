---
type: C++ Function
title: TrimHigher
resource: modeler/internals/skin_intersection.h#L185-L211
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/modeler/internals/skin_intersection/SignedDistanceUnder
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetPositiveSubdivisions
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static std::vector<GeometryPointerType> TrimHigher( const std::vector<GeometryPointerType>& rGeometries, const Point& rLowPoint, std::size_t Component)`

# Calls

- [SignedDistanceUnder](../../../../functions/modeler/internals/skin_intersection/SignedDistanceUnder.md)
- [GetPositiveSubdivisions](../../../../functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetPositiveSubdivisions.md)