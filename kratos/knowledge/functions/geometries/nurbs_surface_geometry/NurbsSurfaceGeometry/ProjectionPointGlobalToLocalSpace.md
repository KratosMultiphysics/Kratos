---
type: C++ Method
title: ProjectionPointGlobalToLocalSpace
resource: geometries/nurbs_surface_geometry.h#L802-L827
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/nurbs_utilities/projection_nurbs_geometry_utilities/NewtonRaphsonSurface
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/nurbs_utilities/projection_nurbs_geometry_utilities/LevenbergMarquardtSurface
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`int ProjectionPointGlobalToLocalSpace( const CoordinatesArrayType& rPointGlobalCoordinates, CoordinatesArrayType& rProjectedPointLocalCoordinates, const double Tolerance = std::numeric_limits<double>::epsilon() ) const override`

# Calls

- [NewtonRaphsonSurface](../../../../functions/utilities/nurbs_utilities/projection_nurbs_geometry_utilities/NewtonRaphsonSurface.md)
- [LevenbergMarquardtSurface](../../../../functions/utilities/nurbs_utilities/projection_nurbs_geometry_utilities/LevenbergMarquardtSurface.md)