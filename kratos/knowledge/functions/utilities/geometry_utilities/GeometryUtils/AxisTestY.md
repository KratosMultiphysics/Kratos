---
type: C++ Method
title: AxisTestY
resource: utilities/geometry_utilities.cpp#L837-L856
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

`bool GeometryUtils::AxisTestY( const double EdgeX, const double EdgeY, const double AbsEdgeX, const double AbsEdgeZ, const array_1d<double,3>& rVertA, const array_1d<double,3>& rVertC, const Point& rBoxHalfSize )`

# Called by

- [TriangleBoxOverlap](../../../../functions/utilities/geometry_utilities/GeometryUtils/TriangleBoxOverlap.md)