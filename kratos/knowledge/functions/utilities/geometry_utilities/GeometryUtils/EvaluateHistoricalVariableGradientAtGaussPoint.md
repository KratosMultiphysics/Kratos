---
type: C++ Method
title: EvaluateHistoricalVariableGradientAtGaussPoint
resource: utilities/geometry_utilities.cpp#L465-L482
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void GeometryUtils::EvaluateHistoricalVariableGradientAtGaussPoint( array_1d<double, 3>& rOutput, const GeometryType& rGeometry, const Variable<double>& rVariable, const Matrix& rGaussPointShapeFunctionDerivativeValues, const int Step)`

# Calls

- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)