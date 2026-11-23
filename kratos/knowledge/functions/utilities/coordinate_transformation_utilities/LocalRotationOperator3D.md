---
type: C++ Function
title: LocalRotationOperator3D
resource: utilities/coordinate_transformation_utilities.h#L959-L1006
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/coordinate_transformation_utilities/Normalize
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void LocalRotationOperator3D( BoundedMatrix<double,TBlockSize,TBlockSize>& rRot, GeometryType::PointType& rThisPoint) const`

# Calls

- [Normalize](../../../functions/utilities/coordinate_transformation_utilities/Normalize.md)