---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/geometries/test_triangle_2d_10.cpp#L70-L73
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/cpp_tests/geometries/test_triangle_2d_10/GenerateReferenceTriangle2D10
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_triangle_2d_10/GenerateCurvedTriangle2D10
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

`KRATOS_TEST_CASE_IN_SUITE(Triangle2D10EdgesNumber, KratosCoreGeometriesFastSuite)`

# Calls

- [GenerateReferenceTriangle2D10](../../../../../functions/tests/cpp_tests/geometries/test_triangle_2d_10/GenerateReferenceTriangle2D10.md)
- [GenerateCurvedTriangle2D10](../../../../../functions/tests/cpp_tests/geometries/test_triangle_2d_10/GenerateCurvedTriangle2D10.md)
- [CrossCheckShapeFunctionsValues](../../../../../functions/tests/cpp_tests/sources/cross_check_shape_functions_values/CrossCheckShapeFunctionsValues.md)
- [TestAllShapeFunctionsLocalGradients](../../../../../functions/tests/cpp_tests/sources/test_shape_function_derivatives/TestAllShapeFunctionsLocalGradients.md)