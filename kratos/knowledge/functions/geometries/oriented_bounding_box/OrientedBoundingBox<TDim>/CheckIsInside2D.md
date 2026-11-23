---
type: C++ Method
title: CheckIsInside2D
resource: geometries/oriented_bounding_box.cpp#L616-L622
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/oriented_bounding_box/OrientedBoundingBox<TDim>/RotateNode2D
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/geometries/oriented_bounding_box/OrientedBoundingBox<2>/IsInside
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool OrientedBoundingBox<TDim>::CheckIsInside2D(array_1d<double, 3>& rCoords) const`

# Calls

- [RotateNode2D](../../../../functions/geometries/oriented_bounding_box/OrientedBoundingBox<TDim>/RotateNode2D.md)

# Called by

- [IsInside](../../../../functions/geometries/oriented_bounding_box/OrientedBoundingBox<2>/IsInside.md)