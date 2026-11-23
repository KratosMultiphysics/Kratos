---
type: C++ Function
title: ApplySlipCondition
resource: utilities/coordinate_transformation_utilities.h#L549-L587
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/coordinate_transformation_utilities/IsSlip
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/coordinate_transformation_utilities/Normalize
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual void ApplySlipCondition(TLocalMatrixType& rLocalMatrix, TLocalVectorType& rLocalVector, GeometryType& rGeometry) const`

# Calls

- [PointsNumber](../../../functions/geometries/geometry/PointsNumber.md)
- [IsSlip](../../../functions/utilities/coordinate_transformation_utilities/IsSlip.md)
- [Normalize](../../../functions/utilities/coordinate_transformation_utilities/Normalize.md)