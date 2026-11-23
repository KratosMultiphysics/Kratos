---
type: C++ Method
title: IsInsideLocalSpace
resource: geometries/nurbs_curve_geometry.h#L313-L333
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/nurbs_curve_geometry/NurbsCurveGeometry/NumberOfKnots
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`int IsInsideLocalSpace( const CoordinatesArrayType& rPointLocalCoordinates, const double Tolerance = std::numeric_limits<double>::epsilon() ) const override`

# Calls

- [NumberOfKnots](../../../../functions/geometries/nurbs_curve_geometry/NurbsCurveGeometry/NumberOfKnots.md)