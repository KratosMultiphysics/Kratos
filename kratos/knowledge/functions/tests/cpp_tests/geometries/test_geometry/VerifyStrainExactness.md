---
type: C++ Function
title: VerifyStrainExactness
resource: tests/cpp_tests/geometries/test_geometry.cpp#L290-L423
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
  - target: functions/includes/amatrix_interface/norm_frobenius
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void VerifyStrainExactness( GeometryType& ThisGeometry, GeometryType::IntegrationMethod ThisMethod )`

# Calls

- [PointsNumber](../../../../../functions/geometries/geometry/PointsNumber.md)
- [prod](../../../../../functions/includes/amatrix_interface/prod.md)
- [norm_frobenius](../../../../../functions/includes/amatrix_interface/norm_frobenius.md)
- [norm_2](../../../../../functions/includes/amatrix_interface/norm_2.md)