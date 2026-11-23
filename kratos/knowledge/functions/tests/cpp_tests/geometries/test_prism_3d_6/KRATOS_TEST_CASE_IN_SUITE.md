---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/geometries/test_prism_3d_6.cpp#L63-L67
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/cpp_tests/geometries/test_prism_3d_6/GenerateRegularPrism3D6
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_geometry/CalculateAreaByIntegration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/CalculateDistanceFrom3DGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(Prism3D6EdgesNumber, KratosCoreGeometriesFastSuite)`

# Calls

- [GenerateRegularPrism3D6](../../../../../functions/tests/cpp_tests/geometries/test_prism_3d_6/GenerateRegularPrism3D6.md)
- [CalculateAreaByIntegration](../../../../../functions/tests/cpp_tests/geometries/test_geometry/CalculateAreaByIntegration.md)
- [CalculateDistanceFrom3DGeometry](../../../../../functions/utilities/geometry_utilities/CalculateDistanceFrom3DGeometry.md)