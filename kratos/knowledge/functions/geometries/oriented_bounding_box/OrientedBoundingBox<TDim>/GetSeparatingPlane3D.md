---
type: C++ Method
title: GetSeparatingPlane3D
resource: geometries/oriented_bounding_box.cpp#L460-L475
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/geometries/oriented_bounding_box/OrientedBoundingBox<3>/SeparatingAxisTheoremHasIntersection
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/oriented_bounding_box/OrientedBoundingBox<TDim>/GetSeparatingPlane
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool OrientedBoundingBox<TDim>::GetSeparatingPlane3D( const array_1d<double, 3>& rRelativePosition, const array_1d<double, 3>& rPlane, const OrientedBoundingBox& rOtherOrientedBoundingBox ) const`

# Called by

- [SeparatingAxisTheoremHasIntersection](../../../../functions/geometries/oriented_bounding_box/OrientedBoundingBox<3>/SeparatingAxisTheoremHasIntersection.md)
- [GetSeparatingPlane](../../../../functions/geometries/oriented_bounding_box/OrientedBoundingBox<TDim>/GetSeparatingPlane.md)