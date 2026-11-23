---
type: C++ Function
title: ProjectionPointGlobalToLocalSpace
resource: geometries/line_2d_2.h#L1156-L1171
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/geometrical_projection_utilities/FastProjectOnLine2D
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`int ProjectionPointGlobalToLocalSpace( const CoordinatesArrayType& rPointGlobalCoordinates, CoordinatesArrayType& rProjectionPointLocalCoordinates, const double Tolerance = std::numeric_limits<double>::epsilon() ) const override`

# Calls

- [FastProjectOnLine2D](../../../functions/utilities/geometrical_projection_utilities/FastProjectOnLine2D.md)