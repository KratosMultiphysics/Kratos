---
type: C++ Method
title: ComputeNurbsShapeFunctionValuesAtSpan
resource: geometries/nurbs_shape_function_utilities/nurbs_surface_shape_functions.h#L298-L377
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_utilities/GetNumberOfControlPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_utilities/GetBinomCoefficient
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ComputeNurbsShapeFunctionValuesAtSpan( const Vector& rKnotsU, const Vector& rKnotsV, const IndexType SpanU, const IndexType SpanV, const Vector& Weights, const double ParameterU, const double ParameterV)`

# Calls

- [GetNumberOfControlPoints](../../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_utilities/GetNumberOfControlPoints.md)
- [GetBinomCoefficient](../../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_utilities/GetBinomCoefficient.md)