---
type: C++ Function
title: NewtonRaphsonSurface
resource: utilities/nurbs_utilities/projection_nurbs_geometry_utilities.h#L284-L397
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/geometries/nurbs_surface_geometry/NurbsSurfaceGeometry/ProjectionPointGlobalToLocalSpace
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_projection_nurbs_geometry/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static bool NewtonRaphsonSurface( CoordinatesArrayType& rProjectedPointLocalCoordinates, const CoordinatesArrayType& rPointGlobalCoordinates, CoordinatesArrayType& rProjectedPointGlobalCoordinates, const NurbsSurfaceGeometry<TDimension, TPointType>& rNurbsSurface, const int MaxIterations = 20, const double Accuracy = 1e-6)`

# Calls

- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)

# Called by

- [ProjectionPointGlobalToLocalSpace](../../../../functions/geometries/nurbs_surface_geometry/NurbsSurfaceGeometry/ProjectionPointGlobalToLocalSpace.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/geometries/test_projection_nurbs_geometry/KRATOS_TEST_CASE_IN_SUITE.md)