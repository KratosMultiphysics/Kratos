---
type: C++ Method
title: ProjectionPointGlobalToLocalSpace
resource: geometries/nurbs_curve_geometry.h#L387-L410
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/nurbs_utilities/projection_nurbs_geometry_utilities/NewtonRaphsonCurve
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/nurbs_utilities/projection_nurbs_geometry_utilities/LevenbergMarquardtCurve
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`int ProjectionPointGlobalToLocalSpace( const CoordinatesArrayType& rPointGlobalCoordinates, CoordinatesArrayType& rProjectedPointLocalCoordinates, const double Tolerance = std::numeric_limits<double>::epsilon() ) const override`

# Calls

- [NewtonRaphsonCurve](../../../../functions/utilities/nurbs_utilities/projection_nurbs_geometry_utilities/NewtonRaphsonCurve.md)
- [LevenbergMarquardtCurve](../../../../functions/utilities/nurbs_utilities/projection_nurbs_geometry_utilities/LevenbergMarquardtCurve.md)