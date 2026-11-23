---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/geometries/test_line_2d_3.cpp#L83-L86
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/cpp_tests/geometries/test_line_2d_3/GeneratePointsUnitXDirectionLine2D3
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_line_2d_3/GeneratePointsDiagonalLine2D3
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_line_2d_3/GeneratePointsParabolaLine2D3
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

`KRATOS_TEST_CASE_IN_SUITE(Line2D3EdgesNumber, KratosCoreGeometriesFastSuite)`

# Calls

- [GeneratePointsUnitXDirectionLine2D3](../../../../../functions/tests/cpp_tests/geometries/test_line_2d_3/GeneratePointsUnitXDirectionLine2D3.md)
- [PointsNumber](../../../../../functions/geometries/geometry/PointsNumber.md)
- [GeneratePointsDiagonalLine2D3](../../../../../functions/tests/cpp_tests/geometries/test_line_2d_3/GeneratePointsDiagonalLine2D3.md)
- [GeneratePointsParabolaLine2D3](../../../../../functions/tests/cpp_tests/geometries/test_line_2d_3/GeneratePointsParabolaLine2D3.md)
- [CrossCheckShapeFunctionsValues](../../../../../functions/tests/cpp_tests/sources/cross_check_shape_functions_values/CrossCheckShapeFunctionsValues.md)
- [TestAllShapeFunctionsLocalGradients](../../../../../functions/tests/cpp_tests/sources/test_shape_function_derivatives/TestAllShapeFunctionsLocalGradients.md)