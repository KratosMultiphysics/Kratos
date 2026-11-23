---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/geometries/test_nurbs_curve_on_surface_geometry.cpp#L367-L395
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/cpp_tests/geometries/test_nurbs_curve_on_surface_geometry/GenerateReferenceBSplineCurveOnBSplineSurface3d
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_nurbs_curve_on_surface_geometry/GenerateReferenceNurbsCurveOnNurbsSurface3d
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_nurbs_curve_on_surface_geometry/GenerateLineNurbsCurveOnSurface3d
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_nurbs_curve_on_surface_geometry/GenerateReferenceNurbsCOS3dforKnotIntersections
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_point/IntegrationPoint/Weight
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/GetDefaultIntegrationMethod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_nurbs_curve_on_surface_geometry/GenerateReferenceNurbs2dforKnotIntersections
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_curve_on_surface_geometry/NurbsCurveOnSurfaceGeometry/CreateQuadraturePointGeometriesSBM
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(BSplineCurveOnSurfaceBSpline, KratosCoreNurbsGeometriesFastSuite)`

# Calls

- [GenerateReferenceBSplineCurveOnBSplineSurface3d](../../../../../functions/tests/cpp_tests/geometries/test_nurbs_curve_on_surface_geometry/GenerateReferenceBSplineCurveOnBSplineSurface3d.md)
- [GenerateReferenceNurbsCurveOnNurbsSurface3d](../../../../../functions/tests/cpp_tests/geometries/test_nurbs_curve_on_surface_geometry/GenerateReferenceNurbsCurveOnNurbsSurface3d.md)
- [GenerateLineNurbsCurveOnSurface3d](../../../../../functions/tests/cpp_tests/geometries/test_nurbs_curve_on_surface_geometry/GenerateLineNurbsCurveOnSurface3d.md)
- [GenerateReferenceNurbsCOS3dforKnotIntersections](../../../../../functions/tests/cpp_tests/geometries/test_nurbs_curve_on_surface_geometry/GenerateReferenceNurbsCOS3dforKnotIntersections.md)
- [Weight](../../../../../functions/integration/integration_point/IntegrationPoint/Weight.md)
- [GetDefaultIntegrationMethod](../../../../../functions/geometries/geometry/GetDefaultIntegrationMethod.md)
- [GenerateReferenceNurbs2dforKnotIntersections](../../../../../functions/tests/cpp_tests/geometries/test_nurbs_curve_on_surface_geometry/GenerateReferenceNurbs2dforKnotIntersections.md)
- [CreateQuadraturePointGeometriesSBM](../../../../../functions/geometries/nurbs_curve_on_surface_geometry/NurbsCurveOnSurfaceGeometry/CreateQuadraturePointGeometriesSBM.md)
- [norm_2](../../../../../functions/includes/amatrix_interface/norm_2.md)