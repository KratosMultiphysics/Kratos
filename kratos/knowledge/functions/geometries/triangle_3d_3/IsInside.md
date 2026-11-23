---
type: C++ Function
title: IsInside
resource: geometries/triangle_3d_3.h#L891-L929
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/geometrical_projection_utilities/FastProject
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool IsInside( const CoordinatesArrayType& rPoint, CoordinatesArrayType& rResult, const double Tolerance = std::numeric_limits<double>::epsilon() ) const override`

# Calls

- [FastProject](../../../functions/utilities/geometrical_projection_utilities/FastProject.md)