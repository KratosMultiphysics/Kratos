---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/geometries/test_quadrilateral_2d_8.cpp#L52-L55
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/cpp_tests/geometries/test_quadrilateral_2d_8/GenerateReferenceQuadrilateral2D8
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

`KRATOS_TEST_CASE_IN_SUITE(Quadrilateral2D8EdgesNumber, KratosCoreGeometriesFastSuite)`

# Calls

- [GenerateReferenceQuadrilateral2D8](../../../../../functions/tests/cpp_tests/geometries/test_quadrilateral_2d_8/GenerateReferenceQuadrilateral2D8.md)
- [CrossCheckShapeFunctionsValues](../../../../../functions/tests/cpp_tests/sources/cross_check_shape_functions_values/CrossCheckShapeFunctionsValues.md)
- [TestAllShapeFunctionsLocalGradients](../../../../../functions/tests/cpp_tests/sources/test_shape_function_derivatives/TestAllShapeFunctionsLocalGradients.md)