---
type: C++ Function
title: GetUpperSpan
resource: geometries/nurbs_shape_function_utilities/nurbs_utilities.h#L49-L57
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_curve_refinement_utilities/NurbsCurveRefinementUtilities/KnotRefinement
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_curve_shape_functions/NurbsCurveShapeFunction/ComputeBSplineShapeFunctionValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_curve_shape_functions/NurbsCurveShapeFunction/ComputeNurbsShapeFunctionValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_surface_refinement_utilities/NurbsSurfaceRefinementUtilities/KnotRefinementU
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_surface_refinement_utilities/NurbsSurfaceRefinementUtilities/KnotRefinementV
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_nurbs_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`inline IndexType GetUpperSpan( const SizeType PolynomialDegree, const Vector& rKnots, const double ParameterT)`

# Called by

- [KnotRefinement](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_curve_refinement_utilities/NurbsCurveRefinementUtilities/KnotRefinement.md)
- [ComputeBSplineShapeFunctionValues](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_curve_shape_functions/NurbsCurveShapeFunction/ComputeBSplineShapeFunctionValues.md)
- [ComputeNurbsShapeFunctionValues](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_curve_shape_functions/NurbsCurveShapeFunction/ComputeNurbsShapeFunctionValues.md)
- [KnotRefinementU](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_surface_refinement_utilities/NurbsSurfaceRefinementUtilities/KnotRefinementU.md)
- [KnotRefinementV](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_surface_refinement_utilities/NurbsSurfaceRefinementUtilities/KnotRefinementV.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/geometries/test_nurbs_utilities/KRATOS_TEST_CASE_IN_SUITE.md)