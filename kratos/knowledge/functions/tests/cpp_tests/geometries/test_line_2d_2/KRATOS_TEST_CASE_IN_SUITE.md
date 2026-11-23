---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/geometries/test_line_2d_2.cpp#L70-L73
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/cpp_tests/geometries/test_line_2d_2/GeneratePointsUnitXDirectionLine2D2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_line_2d_2/GeneratePointsDiagonalLine2D2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometrical_projection_utilities/FastProjectOnLine2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/cross_check_shape_functions_values/CrossCheckShapeFunctionsValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_shape_function_derivatives/TestAllShapeFunctionsLocalGradients
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_line_2d_2/GeneratePointsUnitYDirectionLine2D2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_line_2d_2/GenerateLine2D2WithPoints
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(Line2D2EdgesNumber, KratosCoreGeometriesFastSuite)`

# Calls

- [GeneratePointsUnitXDirectionLine2D2](../../../../../functions/tests/cpp_tests/geometries/test_line_2d_2/GeneratePointsUnitXDirectionLine2D2.md)
- [GeneratePointsDiagonalLine2D2](../../../../../functions/tests/cpp_tests/geometries/test_line_2d_2/GeneratePointsDiagonalLine2D2.md)
- [FastProjectOnLine2D](../../../../../functions/utilities/geometrical_projection_utilities/FastProjectOnLine2D.md)
- [CrossCheckShapeFunctionsValues](../../../../../functions/tests/cpp_tests/sources/cross_check_shape_functions_values/CrossCheckShapeFunctionsValues.md)
- [TestAllShapeFunctionsLocalGradients](../../../../../functions/tests/cpp_tests/sources/test_shape_function_derivatives/TestAllShapeFunctionsLocalGradients.md)
- [GeneratePointsUnitYDirectionLine2D2](../../../../../functions/tests/cpp_tests/geometries/test_line_2d_2/GeneratePointsUnitYDirectionLine2D2.md)
- [GenerateLine2D2WithPoints](../../../../../functions/tests/cpp_tests/geometries/test_line_2d_2/GenerateLine2D2WithPoints.md)