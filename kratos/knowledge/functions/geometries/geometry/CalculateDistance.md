---
type: C++ Function
title: CalculateDistance
resource: geometries/geometry.h#L2937-L2954
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/ClosestPointGlobalToLocalSpace
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual double CalculateDistance( const CoordinatesArrayType& rPointGlobalCoordinates, const double Tolerance = std::numeric_limits<double>::epsilon() ) const`

# Calls

- [ClosestPointGlobalToLocalSpace](../../../functions/geometries/geometry/ClosestPointGlobalToLocalSpace.md)
- [norm_2](../../../functions/includes/amatrix_interface/norm_2.md)