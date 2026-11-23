---
type: C++ Method
title: ComputeNurbsShapeFunctionValuesAtSpan
resource: geometries/nurbs_shape_function_utilities/nurbs_curve_shape_functions.h#L303-L340
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_curve_shape_functions/NurbsCurveShapeFunction/GetFirstNonzeroControlPoint
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_utilities/GetBinomCoefficient
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ComputeNurbsShapeFunctionValuesAtSpan( const Vector& rKnots, const IndexType Span, const Vector& rWeights, const double ParameterT)`

# Calls

- [GetFirstNonzeroControlPoint](../../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_curve_shape_functions/NurbsCurveShapeFunction/GetFirstNonzeroControlPoint.md)
- [GetBinomCoefficient](../../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_utilities/GetBinomCoefficient.md)