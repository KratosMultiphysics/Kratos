---
type: C++ Method
title: RotateNode2D
resource: geometries/oriented_bounding_box.cpp#L566-L585
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/geometries/oriented_bounding_box/OrientedBoundingBox<2>/GetEquivalentRotatedGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/oriented_bounding_box/OrientedBoundingBox<TDim>/CheckIsInside2D
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void OrientedBoundingBox<TDim>::RotateNode2D(array_1d<double, 3>& rCoords) const`

# Called by

- [GetEquivalentRotatedGeometry](../../../../functions/geometries/oriented_bounding_box/OrientedBoundingBox<2>/GetEquivalentRotatedGeometry.md)
- [CheckIsInside2D](../../../../functions/geometries/oriented_bounding_box/OrientedBoundingBox<TDim>/CheckIsInside2D.md)