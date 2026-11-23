---
type: C++ Method
title: GenerateDivision
resource: utilities/divide_triangle_2d_3.cpp#L80-L185
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetNodalDistances
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/split_triangle/TriangleSplit/TriangleSplitMode
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/split_triangle/TriangleSplit/Split_Triangle
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/split_triangle/TriangleSplit/TriangleGetNewConnectivityGID
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void DivideTriangle2D3<TPointType>::GenerateDivision()`

# Calls

- [GetNodalDistances](../../../../functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetNodalDistances.md)
- [TriangleSplitMode](../../../../functions/utilities/split_triangle/TriangleSplit/TriangleSplitMode.md)
- [Split_Triangle](../../../../functions/utilities/split_triangle/TriangleSplit/Split_Triangle.md)
- [TriangleGetNewConnectivityGID](../../../../functions/utilities/split_triangle/TriangleSplit/TriangleGetNewConnectivityGID.md)