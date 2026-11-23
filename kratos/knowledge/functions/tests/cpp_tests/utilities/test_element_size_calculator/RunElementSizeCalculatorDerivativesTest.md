---
type: C++ Function
title: RunElementSizeCalculatorDerivativesTest
resource: tests/cpp_tests/utilities/test_element_size_calculator.cpp#L39-L60
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/utilities/test_element_size_calculator/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void RunElementSizeCalculatorDerivativesTest( Geometry<NodeType>& rGeometry, const TPrimalFunction& rPrimalMethod, const TDerivativeFunction& rDerivativeMethod, const double Delta, const double Tolerance)`

# Calls

- [PointsNumber](../../../../../functions/geometries/geometry/PointsNumber.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/tests/cpp_tests/utilities/test_element_size_calculator/KRATOS_TEST_CASE_IN_SUITE.md)