---
type: C++ Method
title: CreateQuadraturePointGeometriesSBM
resource: geometries/nurbs_curve_on_surface_geometry.h#L691-L828
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/GetDefaultIntegrationMethod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/brep_surface/BrepSurface/KnotsU
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/brep_surface/BrepSurface/KnotsV
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_utilities/GetLowerSpan
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/geometries/brep_curve_on_surface/BrepCurveOnSurface/CreateQuadraturePointGeometries
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_nurbs_curve_on_surface_geometry/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CreateQuadraturePointGeometriesSBM( GeometriesArrayType& rResultGeometries, IndexType NumberOfShapeFunctionDerivatives, const IntegrationPointsArrayType& rIntegrationPoints, IntegrationInfo& rIntegrationInfo)`

# Calls

- [GetDefaultIntegrationMethod](../../../../functions/geometries/geometry/GetDefaultIntegrationMethod.md)
- [KnotsU](../../../../functions/geometries/brep_surface/BrepSurface/KnotsU.md)
- [KnotsV](../../../../functions/geometries/brep_surface/BrepSurface/KnotsV.md)
- [GetLowerSpan](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_utilities/GetLowerSpan.md)

# Called by

- [CreateQuadraturePointGeometries](../../../../functions/geometries/brep_curve_on_surface/BrepCurveOnSurface/CreateQuadraturePointGeometries.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/geometries/test_nurbs_curve_on_surface_geometry/KRATOS_TEST_CASE_IN_SUITE.md)