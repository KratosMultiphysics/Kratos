---
type: C++ Function
title: RotateRHSAux
resource: utilities/coordinate_transformation_utilities.h#L1128-L1238
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
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/coordinate_transformation_utilities/Rotate
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/coordinate_transformation_utilities/RevertRotate
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void RotateRHSAux( TLocalVectorType& rLocalVector, GeometryType& rGeometry, const bool RevertRotation) const`

# Calls

- [PointsNumber](../../../functions/geometries/geometry/PointsNumber.md)
- [IsSlip](../../../functions/utilities/coordinate_transformation_utilities/IsSlip.md)
- [prod](../../../functions/includes/amatrix_interface/prod.md)

# Called by

- [Rotate](../../../functions/utilities/coordinate_transformation_utilities/Rotate.md)
- [RevertRotate](../../../functions/utilities/coordinate_transformation_utilities/RevertRotate.md)