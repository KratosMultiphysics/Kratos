---
type: C++ Method
title: EvaluateHistoricalVariableValueAtGaussPoint
resource: utilities/geometry_utilities.cpp#L416-L437
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/utilities/test_geometry_utils/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void GeometryUtils::EvaluateHistoricalVariableValueAtGaussPoint( TDataType& rOutput, const GeometryType& rGeometry, const Variable<TDataType>& rVariable, const Vector& rGaussPointShapeFunctionValues, const int Step)`

# Calls

- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_geometry_utils/KRATOS_TEST_CASE_IN_SUITE.md)