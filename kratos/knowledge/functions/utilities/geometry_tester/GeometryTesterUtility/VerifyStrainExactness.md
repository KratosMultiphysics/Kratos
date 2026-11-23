---
type: C++ Method
title: VerifyStrainExactness
resource: utilities/geometry_tester.cpp#L735-L857
visibility: private
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

`void GeometryTesterUtility::VerifyStrainExactness( GeometryType& rGeometry, GeometryType::IntegrationMethod ThisMethod, std::stringstream& rErrorMessage )`

# Calls

- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)
- [prod](../../../../functions/includes/amatrix_interface/prod.md)
- [norm_frobenius](../../../../functions/includes/amatrix_interface/norm_frobenius.md)
- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)