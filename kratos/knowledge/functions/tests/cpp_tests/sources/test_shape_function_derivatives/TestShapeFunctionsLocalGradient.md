---
type: C++ Function
title: TestShapeFunctionsLocalGradient
resource: tests/cpp_tests/sources/test_shape_function_derivatives.cpp#L98-L116
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_shape_function_derivatives/FiniteDifference4
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/sources/test_shape_function_derivatives/TestShapeFunctionsLocalGradients
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_shape_function_derivatives/TestShapeFunctionsLocalGradients_IntegrationPointIndex
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_shape_function_derivatives/TestShapeFunctionsLocalGradients_IntegrationPointIndex_ShapeFunctionIndex
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_shape_function_derivatives/TestShapeFunctionsLocalGradients_Point
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void TestShapeFunctionsLocalGradient(Geometry<Node> const& rGeom, Geometry<Node>::IntegrationPointType Point, Matrix const& rLocalGradient)`

# Calls

- [PointsNumber](../../../../../functions/geometries/geometry/PointsNumber.md)
- [FiniteDifference4](../../../../../functions/tests/cpp_tests/sources/test_shape_function_derivatives/FiniteDifference4.md)

# Called by

- [TestShapeFunctionsLocalGradients](../../../../../functions/tests/cpp_tests/sources/test_shape_function_derivatives/TestShapeFunctionsLocalGradients.md)
- [TestShapeFunctionsLocalGradients_IntegrationPointIndex](../../../../../functions/tests/cpp_tests/sources/test_shape_function_derivatives/TestShapeFunctionsLocalGradients_IntegrationPointIndex.md)
- [TestShapeFunctionsLocalGradients_IntegrationPointIndex_ShapeFunctionIndex](../../../../../functions/tests/cpp_tests/sources/test_shape_function_derivatives/TestShapeFunctionsLocalGradients_IntegrationPointIndex_ShapeFunctionIndex.md)
- [TestShapeFunctionsLocalGradients_Point](../../../../../functions/tests/cpp_tests/sources/test_shape_function_derivatives/TestShapeFunctionsLocalGradients_Point.md)