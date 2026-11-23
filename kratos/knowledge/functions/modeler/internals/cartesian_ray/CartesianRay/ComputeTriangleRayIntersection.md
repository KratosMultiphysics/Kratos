---
type: C++ Method
title: ComputeTriangleRayIntersection
resource: modeler/internals/cartesian_ray.h#L309-L346
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/modeler/internals/cartesian_ray/CartesianRay/AddIntersection
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static int ComputeTriangleRayIntersection( const TGeometryType& rTriangleGeometry, const array_1d<double,3>& rLinePoint1, const array_1d<double,3>& rLinePoint2, array_1d<double,3>& rIntersectionPoint, const double epsilon = 1e-12)`

# Called by

- [AddIntersection](../../../../../functions/modeler/internals/cartesian_ray/CartesianRay/AddIntersection.md)