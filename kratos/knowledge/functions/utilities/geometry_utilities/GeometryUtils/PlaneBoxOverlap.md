---
type: C++ Method
title: PlaneBoxOverlap
resource: utilities/geometry_utilities.cpp#L885-L905
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/geometry_utilities/GeometryUtils/TriangleBoxOverlap
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool GeometryUtils::PlaneBoxOverlap( const array_1d<double,3>& rNormal, const double Distance, const array_1d<double,3>& rMaxBox )`

# Called by

- [TriangleBoxOverlap](../../../../functions/utilities/geometry_utilities/GeometryUtils/TriangleBoxOverlap.md)