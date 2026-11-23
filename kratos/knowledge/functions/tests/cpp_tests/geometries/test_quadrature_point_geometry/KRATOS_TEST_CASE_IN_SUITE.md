---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/geometries/test_quadrature_point_geometry.cpp#L112-L133
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/cpp_tests/geometries/test_quadrature_point_geometry/GeneratePointsTriangle2D3
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_quadrature_point_geometry/GenerateQuadraturePointGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_quadrature_point_geometry/GeneratePointsTriangle2D3Update
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/quadrature_points_utility/CreateQuadraturePointsUtility/UpdateFromLocalCoordinates
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_point/IntegrationPoint/Weight
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_quadrature_point_geometry/GenerateQuadraturePointGeometry2
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(QuadraturePointUtility, KratosCoreFastSuite)`

# Calls

- [GeneratePointsTriangle2D3](../../../../../functions/tests/cpp_tests/geometries/test_quadrature_point_geometry/GeneratePointsTriangle2D3.md)
- [GenerateQuadraturePointGeometry](../../../../../functions/tests/cpp_tests/geometries/test_quadrature_point_geometry/GenerateQuadraturePointGeometry.md)
- [GeneratePointsTriangle2D3Update](../../../../../functions/tests/cpp_tests/geometries/test_quadrature_point_geometry/GeneratePointsTriangle2D3Update.md)
- [UpdateFromLocalCoordinates](../../../../../functions/utilities/quadrature_points_utility/CreateQuadraturePointsUtility/UpdateFromLocalCoordinates.md)
- [Weight](../../../../../functions/integration/integration_point/IntegrationPoint/Weight.md)
- [GenerateQuadraturePointGeometry2](../../../../../functions/tests/cpp_tests/geometries/test_quadrature_point_geometry/GenerateQuadraturePointGeometry2.md)