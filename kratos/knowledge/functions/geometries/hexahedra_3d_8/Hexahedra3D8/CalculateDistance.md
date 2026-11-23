---
type: C++ Method
title: CalculateDistance
resource: geometries/hexahedra_3d_8.h#L1115-L1139
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/geometry_utilities/GeometryUtils/PointDistanceToQuadrilateral3D
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`double CalculateDistance( const CoordinatesArrayType& rPointGlobalCoordinates, const double Tolerance = std::numeric_limits<double>::epsilon() ) const override`

# Calls

- [PointDistanceToQuadrilateral3D](../../../../functions/utilities/geometry_utilities/GeometryUtils/PointDistanceToQuadrilateral3D.md)