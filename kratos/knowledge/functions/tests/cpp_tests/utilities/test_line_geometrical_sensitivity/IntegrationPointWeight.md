---
type: C++ Function
title: IntegrationPointWeight
resource: tests/cpp_tests/utilities/test_line_geometrical_sensitivity.cpp#L98-L107
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/cpp_tests/utilities/test_line_geometrical_sensitivity/GetJacobian
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_point/IntegrationPoint/Weight
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/utilities/test_line_geometrical_sensitivity/CheckIntegrationPointWeightSensitivity
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`double IntegrationPointWeight( const Geometry<Point>& rGeometry, GeometryData::IntegrationMethod Quadrature, unsigned int GaussPointIndex)`

# Calls

- [GetJacobian](../../../../../functions/tests/cpp_tests/utilities/test_line_geometrical_sensitivity/GetJacobian.md)
- [prod](../../../../../functions/includes/amatrix_interface/prod.md)
- [Weight](../../../../../functions/integration/integration_point/IntegrationPoint/Weight.md)

# Called by

- [CheckIntegrationPointWeightSensitivity](../../../../../functions/tests/cpp_tests/utilities/test_line_geometrical_sensitivity/CheckIntegrationPointWeightSensitivity.md)