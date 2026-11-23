---
type: C++ Function
title: Triangle2D3AusasIncisedShapeFunctions
resource: modified_shape_functions/triangle_2d_3_ausas_incised_shape_functions.h#L209-L212
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetNodalDistances
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`Triangle2D3AusasIncisedShapeFunctions(Triangle2D3AusasIncisedShapeFunctions const& rOther) : Triangle2D3AusasModifiedShapeFunctions(rOther.GetInputGeometry(), rOther.GetNodalDistances()), mExtraEdgeRatios(rOther.mExtraEdgeRatios)`

# Calls

- [GetNodalDistances](../../../functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetNodalDistances.md)