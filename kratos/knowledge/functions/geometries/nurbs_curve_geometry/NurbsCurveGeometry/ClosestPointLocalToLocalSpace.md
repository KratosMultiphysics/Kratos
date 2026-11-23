---
type: C++ Method
title: ClosestPointLocalToLocalSpace
resource: geometries/nurbs_curve_geometry.h#L347-L373
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/nurbs_curve_geometry/NurbsCurveGeometry/NumberOfKnots
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual int ClosestPointLocalToLocalSpace( const CoordinatesArrayType& rPointLocalCoordinates, CoordinatesArrayType& rClosestPointLocalCoordinates, const double Tolerance = std::numeric_limits<double>::epsilon() ) const override`

# Calls

- [NumberOfKnots](../../../../functions/geometries/nurbs_curve_geometry/NurbsCurveGeometry/NumberOfKnots.md)