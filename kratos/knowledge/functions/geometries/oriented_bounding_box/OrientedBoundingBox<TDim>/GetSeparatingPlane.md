---
type: C++ Method
title: GetSeparatingPlane
resource: geometries/oriented_bounding_box.cpp#L423-L435
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/oriented_bounding_box/OrientedBoundingBox<TDim>/GetSeparatingPlane2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/oriented_bounding_box/OrientedBoundingBox<TDim>/GetSeparatingPlane3D
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool OrientedBoundingBox<TDim>::GetSeparatingPlane( const array_1d<double, 3>& rRelativePosition, const array_1d<double, 3>& rPlane, const OrientedBoundingBox& rOtherOrientedBoundingBox ) const`

# Calls

- [GetSeparatingPlane2D](../../../../functions/geometries/oriented_bounding_box/OrientedBoundingBox<TDim>/GetSeparatingPlane2D.md)
- [GetSeparatingPlane3D](../../../../functions/geometries/oriented_bounding_box/OrientedBoundingBox<TDim>/GetSeparatingPlane3D.md)