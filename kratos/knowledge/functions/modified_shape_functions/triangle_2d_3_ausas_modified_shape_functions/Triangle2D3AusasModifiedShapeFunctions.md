---
type: C++ Function
title: Triangle2D3AusasModifiedShapeFunctions
resource: modified_shape_functions/triangle_2d_3_ausas_modified_shape_functions.h#L182-L189
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetNodalDistances
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`Triangle2D3AusasModifiedShapeFunctions(Triangle2D3AusasModifiedShapeFunctions const& rOther) : AusasModifiedShapeFunctions(rOther.GetInputGeometry(), rOther.GetNodalDistances()), mpTriangleSplitter(new DivideTriangle2D3<Node>(*rOther.GetInputGeometry(), rOther.GetNodalDistances()))`

# Calls

- [GetNodalDistances](../../../functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetNodalDistances.md)