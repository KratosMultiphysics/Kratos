---
type: C++ Method
title: CheckIsInside3D
resource: geometries/oriented_bounding_box.cpp#L628-L637
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/oriented_bounding_box/OrientedBoundingBox<TDim>/RotateNode3D
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/geometries/oriented_bounding_box/OrientedBoundingBox<3>/IsInside
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool OrientedBoundingBox<TDim>::CheckIsInside3D( array_1d<double, 3>& rCoords, BoundedMatrix<double, 4, 4> rInvertedRotationMatrix ) const`

# Calls

- [RotateNode3D](../../../../functions/geometries/oriented_bounding_box/OrientedBoundingBox<TDim>/RotateNode3D.md)

# Called by

- [IsInside](../../../../functions/geometries/oriented_bounding_box/OrientedBoundingBox<3>/IsInside.md)