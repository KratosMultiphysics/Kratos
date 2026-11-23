---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/geometries/test_quadrilateral_3d_9.cpp#L66-L71
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/cpp_tests/geometries/test_quadrilateral_3d_9/CheckSamePoint
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/cross_check_shape_functions_values/CrossCheckShapeFunctionsValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_shape_function_derivatives/TestAllShapeFunctionsLocalGradients
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/math_utils/Norm3
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D9EdgesNumber, KratosCoreGeometriesFastSuite)`

# Calls

- [CheckSamePoint](../../../../../functions/tests/cpp_tests/geometries/test_quadrilateral_3d_9/CheckSamePoint.md)
- [CrossCheckShapeFunctionsValues](../../../../../functions/tests/cpp_tests/sources/cross_check_shape_functions_values/CrossCheckShapeFunctionsValues.md)
- [TestAllShapeFunctionsLocalGradients](../../../../../functions/tests/cpp_tests/sources/test_shape_function_derivatives/TestAllShapeFunctionsLocalGradients.md)
- [Norm3](../../../../../functions/utilities/math_utils/Norm3.md)