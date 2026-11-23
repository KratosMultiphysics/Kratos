---
type: C++ Function
title: CheckIntegrationPointWeightSensitivity
resource: tests/cpp_tests/utilities/test_line_geometrical_sensitivity.cpp#L109-L145
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/cpp_tests/utilities/test_line_geometrical_sensitivity/GetJacobian
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_line_geometrical_sensitivity/IntegrationPointWeight
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_point/IntegrationPoint/Weight
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/utilities/test_line_geometrical_sensitivity/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CheckIntegrationPointWeightSensitivity( Geometry<Point>::Pointer pGeometry, GeometryData::IntegrationMethod Quadrature, double Perturbation = 1e-7, double Tolerance = 1e-7)`

# Calls

- [GetJacobian](../../../../../functions/tests/cpp_tests/utilities/test_line_geometrical_sensitivity/GetJacobian.md)
- [IntegrationPointWeight](../../../../../functions/tests/cpp_tests/utilities/test_line_geometrical_sensitivity/IntegrationPointWeight.md)
- [PointsNumber](../../../../../functions/geometries/geometry/PointsNumber.md)
- [Weight](../../../../../functions/integration/integration_point/IntegrationPoint/Weight.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/tests/cpp_tests/utilities/test_line_geometrical_sensitivity/KRATOS_TEST_CASE_IN_SUITE.md)