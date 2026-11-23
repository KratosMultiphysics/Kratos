---
type: C++ Function
title: ClosestPointLocalCoordinates
resource: geometries/geometry.h#L2842-L2849
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/ClosestPointGlobalToLocalSpace
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual int ClosestPointLocalCoordinates( const CoordinatesArrayType& rPointGlobalCoordinates, CoordinatesArrayType& rClosestPointLocalCoordinates, const double Tolerance = std::numeric_limits<double>::epsilon() ) const`

# Calls

- [ClosestPointGlobalToLocalSpace](../../../functions/geometries/geometry/ClosestPointGlobalToLocalSpace.md)