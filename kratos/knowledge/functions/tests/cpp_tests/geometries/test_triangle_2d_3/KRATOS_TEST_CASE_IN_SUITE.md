---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/geometries/test_triangle_2d_3.cpp#L102-L106
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/cpp_tests/geometries/test_triangle_2d_3/GeneratePointsRightTriangle2D3
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_triangle_2d_3/GenerateNodesRightTriangle2D3
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_triangle_2d_3/GeneratePointsIrregularTriangle2D3
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_geometry/CalculateAreaByIntegration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/cross_check_shape_functions_values/CrossCheckShapeFunctionsValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_shape_function_derivatives/TestAllShapeFunctionsLocalGradients
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(Triangle2D3EdgesNumber, KratosCoreGeometriesFastSuite)`

# Calls

- [GeneratePointsRightTriangle2D3](../../../../../functions/tests/cpp_tests/geometries/test_triangle_2d_3/GeneratePointsRightTriangle2D3.md)
- [GenerateNodesRightTriangle2D3](../../../../../functions/tests/cpp_tests/geometries/test_triangle_2d_3/GenerateNodesRightTriangle2D3.md)
- [GeneratePointsIrregularTriangle2D3](../../../../../functions/tests/cpp_tests/geometries/test_triangle_2d_3/GeneratePointsIrregularTriangle2D3.md)
- [CalculateAreaByIntegration](../../../../../functions/tests/cpp_tests/geometries/test_geometry/CalculateAreaByIntegration.md)
- [CrossCheckShapeFunctionsValues](../../../../../functions/tests/cpp_tests/sources/cross_check_shape_functions_values/CrossCheckShapeFunctionsValues.md)
- [TestAllShapeFunctionsLocalGradients](../../../../../functions/tests/cpp_tests/sources/test_shape_function_derivatives/TestAllShapeFunctionsLocalGradients.md)
- [norm_2](../../../../../functions/includes/amatrix_interface/norm_2.md)