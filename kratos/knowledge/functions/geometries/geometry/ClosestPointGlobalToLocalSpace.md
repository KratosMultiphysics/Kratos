---
type: C++ Function
title: ClosestPointGlobalToLocalSpace
resource: geometries/geometry.h#L2899-L2923
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/geometries/geometry/ClosestPoint
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/ClosestPointLocalCoordinates
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/CalculateDistance
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual int ClosestPointGlobalToLocalSpace( const CoordinatesArrayType& rPointGlobalCoordinates, CoordinatesArrayType& rClosestPointLocalCoordinates, const double Tolerance = std::numeric_limits<double>::epsilon() ) const`

# Called by

- [ClosestPoint](../../../functions/geometries/geometry/ClosestPoint.md)
- [ClosestPointLocalCoordinates](../../../functions/geometries/geometry/ClosestPointLocalCoordinates.md)
- [CalculateDistance](../../../functions/geometries/geometry/CalculateDistance.md)