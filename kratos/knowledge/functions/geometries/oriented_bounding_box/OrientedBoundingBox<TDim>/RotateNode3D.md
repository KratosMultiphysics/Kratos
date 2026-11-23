---
type: C++ Method
title: RotateNode3D
resource: geometries/oriented_bounding_box.cpp#L591-L610
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/geometries/oriented_bounding_box/OrientedBoundingBox<3>/GetEquivalentRotatedGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/oriented_bounding_box/OrientedBoundingBox<TDim>/CheckIsInside3D
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void OrientedBoundingBox<TDim>::RotateNode3D( array_1d<double, 3>& rCoords, const BoundedMatrix<double, 4, 4>& rInvertedRotationMatrix ) const`

# Calls

- [prod](../../../../functions/includes/amatrix_interface/prod.md)

# Called by

- [GetEquivalentRotatedGeometry](../../../../functions/geometries/oriented_bounding_box/OrientedBoundingBox<3>/GetEquivalentRotatedGeometry.md)
- [CheckIsInside3D](../../../../functions/geometries/oriented_bounding_box/OrientedBoundingBox<TDim>/CheckIsInside3D.md)