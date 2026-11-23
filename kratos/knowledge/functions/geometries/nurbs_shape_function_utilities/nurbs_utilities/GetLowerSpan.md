---
type: C++ Function
title: GetLowerSpan
resource: geometries/nurbs_shape_function_utilities/nurbs_utilities.h#L63-L81
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/geometries/nurbs_curve_on_surface_geometry/NurbsCurveOnSurfaceGeometry/CreateQuadraturePointGeometriesSBM
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_surface_shape_functions/NurbsSurfaceShapeFunction/ComputeBSplineShapeFunctionValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_surface_shape_functions/NurbsSurfaceShapeFunction/ComputeNurbsShapeFunctionValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_volume_refinement_utilities/NurbsVolumeRefinementUtilities/KnotRefinementU
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_volume_refinement_utilities/NurbsVolumeRefinementUtilities/KnotRefinementV
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_volume_refinement_utilities/NurbsVolumeRefinementUtilities/KnotRefinementW
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_volume_shape_functions/NurbsVolumeShapeFunction/ComputeBSplineShapeFunctionValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_volume_shape_functions/NurbsVolumeShapeFunction/ComputeNurbsShapeFunctionValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_surface_geometry/NurbsSurfaceGeometry/CalculateEstimatedKnotLengthness
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_surface_on_volume_geometry/NurbsSurfaceOnVolumeGeometry/CreateQuadraturePointGeometriesSbm
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_nurbs_curve_shape_functions/GenerateReferenceNurbsCurveShapeFunction
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_nurbs_surface_shape_functions/GenerateReferenceHemisphereSurfaceShapeFunctions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_nurbs_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_nurbs_volume_shape_functions/GenerateVolumeShapeFunctions
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`inline IndexType GetLowerSpan( const SizeType PolynomialDegree, const Vector& rKnots, const double ParameterT, const double tolerance = 1e-12)`

# Called by

- [CreateQuadraturePointGeometriesSBM](../../../../functions/geometries/nurbs_curve_on_surface_geometry/NurbsCurveOnSurfaceGeometry/CreateQuadraturePointGeometriesSBM.md)
- [ComputeBSplineShapeFunctionValues](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_surface_shape_functions/NurbsSurfaceShapeFunction/ComputeBSplineShapeFunctionValues.md)
- [ComputeNurbsShapeFunctionValues](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_surface_shape_functions/NurbsSurfaceShapeFunction/ComputeNurbsShapeFunctionValues.md)
- [KnotRefinementU](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_volume_refinement_utilities/NurbsVolumeRefinementUtilities/KnotRefinementU.md)
- [KnotRefinementV](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_volume_refinement_utilities/NurbsVolumeRefinementUtilities/KnotRefinementV.md)
- [KnotRefinementW](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_volume_refinement_utilities/NurbsVolumeRefinementUtilities/KnotRefinementW.md)
- [ComputeBSplineShapeFunctionValues](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_volume_shape_functions/NurbsVolumeShapeFunction/ComputeBSplineShapeFunctionValues.md)
- [ComputeNurbsShapeFunctionValues](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_volume_shape_functions/NurbsVolumeShapeFunction/ComputeNurbsShapeFunctionValues.md)
- [CalculateEstimatedKnotLengthness](../../../../functions/geometries/nurbs_surface_geometry/NurbsSurfaceGeometry/CalculateEstimatedKnotLengthness.md)
- [CreateQuadraturePointGeometriesSbm](../../../../functions/geometries/nurbs_surface_on_volume_geometry/NurbsSurfaceOnVolumeGeometry/CreateQuadraturePointGeometriesSbm.md)
- [GenerateReferenceNurbsCurveShapeFunction](../../../../functions/tests/cpp_tests/geometries/test_nurbs_curve_shape_functions/GenerateReferenceNurbsCurveShapeFunction.md)
- [GenerateReferenceHemisphereSurfaceShapeFunctions](../../../../functions/tests/cpp_tests/geometries/test_nurbs_surface_shape_functions/GenerateReferenceHemisphereSurfaceShapeFunctions.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/geometries/test_nurbs_utilities/KRATOS_TEST_CASE_IN_SUITE.md)
- [GenerateVolumeShapeFunctions](../../../../functions/tests/cpp_tests/geometries/test_nurbs_volume_shape_functions/GenerateVolumeShapeFunctions.md)