---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/geometries/test_line_2d_5.cpp#L93-L96
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/cpp_tests/geometries/test_line_2d_5/GeneratePointsUnitXDirectionLine2D5
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_line_2d_5/GeneratePointsDiagonalLine2D5
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_line_2d_5/GeneratePointsParabolaLine2D5
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/cross_check_shape_functions_values/CrossCheckShapeFunctionsValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_shape_function_derivatives/TestAllShapeFunctionsLocalGradients
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(Line2D5EdgesNumber, KratosCoreGeometriesFastSuite)`

# Calls

- [GeneratePointsUnitXDirectionLine2D5](../../../../../functions/tests/cpp_tests/geometries/test_line_2d_5/GeneratePointsUnitXDirectionLine2D5.md)
- [PointsNumber](../../../../../functions/geometries/geometry/PointsNumber.md)
- [GeneratePointsDiagonalLine2D5](../../../../../functions/tests/cpp_tests/geometries/test_line_2d_5/GeneratePointsDiagonalLine2D5.md)
- [GeneratePointsParabolaLine2D5](../../../../../functions/tests/cpp_tests/geometries/test_line_2d_5/GeneratePointsParabolaLine2D5.md)
- [CrossCheckShapeFunctionsValues](../../../../../functions/tests/cpp_tests/sources/cross_check_shape_functions_values/CrossCheckShapeFunctionsValues.md)
- [TestAllShapeFunctionsLocalGradients](../../../../../functions/tests/cpp_tests/sources/test_shape_function_derivatives/TestAllShapeFunctionsLocalGradients.md)