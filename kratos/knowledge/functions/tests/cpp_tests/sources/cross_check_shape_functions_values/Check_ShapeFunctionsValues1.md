---
type: C++ Function
title: Check_ShapeFunctionsValues1
resource: tests/cpp_tests/sources/cross_check_shape_functions_values.cpp#L34-L43
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/sources/cross_check_shape_functions_values/CrossCheckShapeFunctionsValues
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Check_ShapeFunctionsValues1(Geometry<Node> const& rGeom)`

# Calls

- [PointsNumber](../../../../../functions/geometries/geometry/PointsNumber.md)

# Called by

- [CrossCheckShapeFunctionsValues](../../../../../functions/tests/cpp_tests/sources/cross_check_shape_functions_values/CrossCheckShapeFunctionsValues.md)