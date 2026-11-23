---
type: C++ Function
title: TestShapeFunctionsLocalGradients_IntegrationPointIndex_ShapeFunctionIndex
resource: tests/cpp_tests/sources/test_shape_function_derivatives.cpp#L54-L65
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_shape_function_derivatives/TestShapeFunctionsLocalGradient
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/sources/test_shape_function_derivatives/TestAllShapeFunctionsLocalGradients
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void TestShapeFunctionsLocalGradients_IntegrationPointIndex_ShapeFunctionIndex( Geometry<Node> const& rGeom, GeometryData::IntegrationMethod ThisMethod)`

# Calls

- [PointsNumber](../../../../../functions/geometries/geometry/PointsNumber.md)
- [TestShapeFunctionsLocalGradient](../../../../../functions/tests/cpp_tests/sources/test_shape_function_derivatives/TestShapeFunctionsLocalGradient.md)

# Called by

- [TestAllShapeFunctionsLocalGradients](../../../../../functions/tests/cpp_tests/sources/test_shape_function_derivatives/TestAllShapeFunctionsLocalGradients.md)