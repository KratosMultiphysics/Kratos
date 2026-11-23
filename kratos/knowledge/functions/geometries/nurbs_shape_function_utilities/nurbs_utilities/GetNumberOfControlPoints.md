---
type: C++ Function
title: GetNumberOfControlPoints
resource: geometries/nurbs_shape_function_utilities/nurbs_utilities.h#L108-L111
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_surface_shape_functions/NurbsSurfaceShapeFunction/ComputeNurbsShapeFunctionValuesAtSpan
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_surface_shape_functions/NurbsSurfaceShapeFunction/GetControlPointIndex
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_volume_shape_functions/NurbsVolumeShapeFunction/GetControlPointIndex
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_surface_geometry/NurbsSurfaceGeometry/CheckAndFitKnotVectors
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/CheckAndFitKnotVectors
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_nurbs_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`inline SizeType GetNumberOfControlPoints(const SizeType PolynomialDegree, const SizeType NumberOfKnots)`

# Called by

- [ComputeNurbsShapeFunctionValuesAtSpan](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_surface_shape_functions/NurbsSurfaceShapeFunction/ComputeNurbsShapeFunctionValuesAtSpan.md)
- [GetControlPointIndex](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_surface_shape_functions/NurbsSurfaceShapeFunction/GetControlPointIndex.md)
- [GetControlPointIndex](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_volume_shape_functions/NurbsVolumeShapeFunction/GetControlPointIndex.md)
- [CheckAndFitKnotVectors](../../../../functions/geometries/nurbs_surface_geometry/NurbsSurfaceGeometry/CheckAndFitKnotVectors.md)
- [CheckAndFitKnotVectors](../../../../functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/CheckAndFitKnotVectors.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/geometries/test_nurbs_utilities/KRATOS_TEST_CASE_IN_SUITE.md)