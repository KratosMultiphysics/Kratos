---
type: C++ Method
title: GenerateDivision
resource: utilities/divide_tetrahedra_3d_4.cpp#L81-L186
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetNodalDistances
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/split_tetrahedra/TetrahedraSplit/TetrahedraSplitMode
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/split_tetrahedra/TetrahedraSplit/Split_Tetrahedra
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/split_tetrahedra/TetrahedraSplit/TetrahedraGetNewConnectivityGID
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void DivideTetrahedra3D4<TPointType>::GenerateDivision()`

# Calls

- [GetNodalDistances](../../../../functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetNodalDistances.md)
- [TetrahedraSplitMode](../../../../functions/utilities/split_tetrahedra/TetrahedraSplit/TetrahedraSplitMode.md)
- [Split_Tetrahedra](../../../../functions/utilities/split_tetrahedra/TetrahedraSplit/Split_Tetrahedra.md)
- [TetrahedraGetNewConnectivityGID](../../../../functions/utilities/split_tetrahedra/TetrahedraSplit/TetrahedraGetNewConnectivityGID.md)