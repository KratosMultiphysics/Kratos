---
type: C++ Function
title: ProjectionPointGlobalToLocalSpace
resource: geometries/quadrilateral_3d_4.h#L1463-L1506
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`int ProjectionPointGlobalToLocalSpace( const CoordinatesArrayType& rPointGlobalCoordinates, CoordinatesArrayType& rProjectionPointLocalCoordinates, const double Tolerance = std::numeric_limits<double>::epsilon() ) const override`

# Calls

- [norm_2](../../../functions/includes/amatrix_interface/norm_2.md)