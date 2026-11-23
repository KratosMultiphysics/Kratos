---
type: C++ Function
title: CalculateTetrahedraIntersectionPoints
resource: utilities/geometry_utilities.h#L556-L591
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/geometry_utilities/CalculateTetrahedraDistances
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/CalculateTriangleDistances
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/CalculateExactDistancesToPlane
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static int CalculateTetrahedraIntersectionPoints( const GeometryType& rGeometry, array_1d<double, TSize1>& rDistances, array_1d<Point, TSize2>& rIntersectionPoints )`

# Called by

- [CalculateTetrahedraDistances](../../../functions/utilities/geometry_utilities/CalculateTetrahedraDistances.md)
- [CalculateTriangleDistances](../../../functions/utilities/geometry_utilities/CalculateTriangleDistances.md)
- [CalculateExactDistancesToPlane](../../../functions/utilities/geometry_utilities/CalculateExactDistancesToPlane.md)