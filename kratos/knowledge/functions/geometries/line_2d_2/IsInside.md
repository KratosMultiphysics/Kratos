---
type: C++ Function
title: IsInside
resource: geometries/line_2d_2.h#L952-L978
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/geometrical_projection_utilities/FastProjectOnLine2D
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool IsInside( const CoordinatesArrayType& rPoint, CoordinatesArrayType& rResult, const double Tolerance = std::numeric_limits<double>::epsilon() ) const override`

# Calls

- [FastProjectOnLine2D](../../../functions/utilities/geometrical_projection_utilities/FastProjectOnLine2D.md)