---
type: C++ Function
title: GetJacobian
resource: tests/cpp_tests/utilities/test_line_geometrical_sensitivity.cpp#L76-L96
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/utilities/test_line_geometrical_sensitivity/IntegrationPointWeight
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_line_geometrical_sensitivity/CheckIntegrationPointWeightSensitivity
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`DenseMatrix<double> GetJacobian( const Geometry<Point>& rGeometry, GeometryData::IntegrationMethod Quadrature, unsigned int GaussPointIndex)`

# Calls

- [PointsNumber](../../../../../functions/geometries/geometry/PointsNumber.md)
- [prod](../../../../../functions/includes/amatrix_interface/prod.md)

# Called by

- [IntegrationPointWeight](../../../../../functions/tests/cpp_tests/utilities/test_line_geometrical_sensitivity/IntegrationPointWeight.md)
- [CheckIntegrationPointWeightSensitivity](../../../../../functions/tests/cpp_tests/utilities/test_line_geometrical_sensitivity/CheckIntegrationPointWeightSensitivity.md)