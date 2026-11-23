---
type: C++ Method
title: TriangleGetNewConnectivityGID
resource: utilities/split_triangle.h#L96-L104
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/divide_triangle_2d_3/DivideTriangle2D3<TPointType>/GenerateDivision
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`inline static void TriangleGetNewConnectivityGID(const int triangle_index, const int t[12], const int aux_ids[6], int* id0, int* id1, int* id2)`

# Called by

- [GenerateDivision](../../../../functions/utilities/divide_triangle_2d_3/DivideTriangle2D3<TPointType>/GenerateDivision.md)