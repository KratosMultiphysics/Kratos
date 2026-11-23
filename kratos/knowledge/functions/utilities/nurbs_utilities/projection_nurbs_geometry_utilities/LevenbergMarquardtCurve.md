---
type: C++ Function
title: LevenbergMarquardtCurve
resource: utilities/nurbs_utilities/projection_nurbs_geometry_utilities.h#L149-L265
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/geometries/brep_curve_on_surface/BrepCurveOnSurface/ProjectionPointGlobalToLocalSpace
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_curve_geometry/NurbsCurveGeometry/ProjectionPointGlobalToLocalSpace
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_curve_on_surface_geometry/NurbsCurveOnSurfaceGeometry/ProjectionPointGlobalToLocalSpace
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_projection_nurbs_geometry/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static bool LevenbergMarquardtCurve( CoordinatesArrayType& rProjectedPointLocalCoordinates, const CoordinatesArrayType& rPointGlobalCoordinates, CoordinatesArrayType& rProjectedPointGlobalCoordinates, const Geometry<TPointType>& rGeometry, const int MaxIterations = 30, const double Accuracy = 1e-6)`

# Calls

- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)

# Called by

- [ProjectionPointGlobalToLocalSpace](../../../../functions/geometries/brep_curve_on_surface/BrepCurveOnSurface/ProjectionPointGlobalToLocalSpace.md)
- [ProjectionPointGlobalToLocalSpace](../../../../functions/geometries/nurbs_curve_geometry/NurbsCurveGeometry/ProjectionPointGlobalToLocalSpace.md)
- [ProjectionPointGlobalToLocalSpace](../../../../functions/geometries/nurbs_curve_on_surface_geometry/NurbsCurveOnSurfaceGeometry/ProjectionPointGlobalToLocalSpace.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/geometries/test_projection_nurbs_geometry/KRATOS_TEST_CASE_IN_SUITE.md)