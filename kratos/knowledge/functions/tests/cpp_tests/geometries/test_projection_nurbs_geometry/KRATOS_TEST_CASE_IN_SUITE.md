---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/geometries/test_projection_nurbs_geometry.cpp#L210-L264
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/cpp_tests/geometries/test_projection_nurbs_geometry/GenerateReferenceCurveForProjection2d
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/nurbs_utilities/projection_nurbs_geometry_utilities/NewtonRaphsonCurve
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_interval/NurbsInterval/MinParameter
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/nurbs_utilities/projection_nurbs_geometry_utilities/LevenbergMarquardtCurve
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_projection_nurbs_geometry/GenerateReferenceUnrefinedQuarterSphereGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/nurbs_utilities/projection_nurbs_geometry_utilities/NewtonRaphsonSurface
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/nurbs_utilities/projection_nurbs_geometry_utilities/LevenbergMarquardtSurface
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_projection_nurbs_geometry/GenerateReferenceC0HemisphereGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(NurbsCurveGeometryProjection2d, KratosCoreNurbsGeometriesFastSuite)`

# Calls

- [GenerateReferenceCurveForProjection2d](../../../../../functions/tests/cpp_tests/geometries/test_projection_nurbs_geometry/GenerateReferenceCurveForProjection2d.md)
- [NewtonRaphsonCurve](../../../../../functions/utilities/nurbs_utilities/projection_nurbs_geometry_utilities/NewtonRaphsonCurve.md)
- [MinParameter](../../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_interval/NurbsInterval/MinParameter.md)
- [LevenbergMarquardtCurve](../../../../../functions/utilities/nurbs_utilities/projection_nurbs_geometry_utilities/LevenbergMarquardtCurve.md)
- [GenerateReferenceUnrefinedQuarterSphereGeometry](../../../../../functions/tests/cpp_tests/geometries/test_projection_nurbs_geometry/GenerateReferenceUnrefinedQuarterSphereGeometry.md)
- [NewtonRaphsonSurface](../../../../../functions/utilities/nurbs_utilities/projection_nurbs_geometry_utilities/NewtonRaphsonSurface.md)
- [LevenbergMarquardtSurface](../../../../../functions/utilities/nurbs_utilities/projection_nurbs_geometry_utilities/LevenbergMarquardtSurface.md)
- [GenerateReferenceC0HemisphereGeometry](../../../../../functions/tests/cpp_tests/geometries/test_projection_nurbs_geometry/GenerateReferenceC0HemisphereGeometry.md)