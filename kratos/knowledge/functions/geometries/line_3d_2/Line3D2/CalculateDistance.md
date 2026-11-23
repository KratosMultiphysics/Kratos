---
type: C++ Method
title: CalculateDistance
resource: geometries/line_3d_2.h#L462-L469
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/geometry_utilities/GeometryUtils/PointDistanceToLineSegment3D
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`double CalculateDistance( const CoordinatesArrayType& rPointGlobalCoordinates, const double Tolerance = std::numeric_limits<double>::epsilon() ) const override`

# Calls

- [PointDistanceToLineSegment3D](../../../../functions/utilities/geometry_utilities/GeometryUtils/PointDistanceToLineSegment3D.md)