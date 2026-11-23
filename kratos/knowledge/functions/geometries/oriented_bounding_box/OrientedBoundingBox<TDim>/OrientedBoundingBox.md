---
type: C++ Method
title: OrientedBoundingBox
resource: geometries/oriented_bounding_box.cpp#L26-L33
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`OrientedBoundingBox<TDim>::OrientedBoundingBox( const array_1d<double, 3>& rCenterCoords, const array_1d<array_1d<double, 3>, TDim>& rOrientationVectors, const array_1d<double, TDim>& rHalfLength ) : mPointCenter(rCenterCoords), mOrientationVectors(rOrientationVectors), mHalfLength(rHalfLength)`

# Calls

- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)