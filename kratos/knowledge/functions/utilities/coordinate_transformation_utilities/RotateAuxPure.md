---
type: C++ Function
title: RotateAuxPure
resource: utilities/coordinate_transformation_utilities.h#L854-L937
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/coordinate_transformation_utilities/IsSlip
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void RotateAuxPure(TLocalMatrixType& rLocalMatrix, TLocalVectorType& rLocalVector, GeometryType& rGeometry) const`

# Calls

- [IsSlip](../../../functions/utilities/coordinate_transformation_utilities/IsSlip.md)
- [prod](../../../functions/includes/amatrix_interface/prod.md)