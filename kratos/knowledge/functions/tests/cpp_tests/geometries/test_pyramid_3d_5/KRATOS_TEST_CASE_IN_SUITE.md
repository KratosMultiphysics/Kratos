---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/geometries/test_pyramid_3d_5.cpp#L43-L49
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/cpp_tests/geometries/test_pyramid_3d_5/GenerateRegularPyramid3D5
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

`KRATOS_TEST_CASE_IN_SUITE(Pyramid3D5GetGeometryType, KratosCoreGeometriesFastSuite)`

# Calls

- [GenerateRegularPyramid3D5](../../../../../functions/tests/cpp_tests/geometries/test_pyramid_3d_5/GenerateRegularPyramid3D5.md)
- [CalculateAreaByIntegration](../../../../../functions/tests/cpp_tests/geometries/test_geometry/CalculateAreaByIntegration.md)
- [CalculateDistanceFrom3DGeometry](../../../../../functions/utilities/geometry_utilities/CalculateDistanceFrom3DGeometry.md)