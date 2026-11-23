---
type: C++ Function
title: Triangle2D3ModifiedShapeFunctions
resource: modified_shape_functions/triangle_2d_3_modified_shape_functions.h#L185-L192
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetNodalDistances
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`Triangle2D3ModifiedShapeFunctions(Triangle2D3ModifiedShapeFunctions const& rOther) : ModifiedShapeFunctions(rOther.GetInputGeometry(), rOther.GetNodalDistances()), mpTriangleSplitter(new DivideTriangle2D3<Node>(*rOther.GetInputGeometry(), rOther.GetNodalDistances()))`

# Calls

- [GetNodalDistances](../../../functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetNodalDistances.md)