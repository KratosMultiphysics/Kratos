---
type: C++ Function
title: Tetrahedra3D4ModifiedShapeFunctions
resource: modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h#L185-L192
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetNodalDistances
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`Tetrahedra3D4ModifiedShapeFunctions(Tetrahedra3D4ModifiedShapeFunctions const& rOther) : ModifiedShapeFunctions(rOther.GetInputGeometry(), rOther.GetNodalDistances()), mpTetrahedraSplitter(new DivideTetrahedra3D4<Node>(*rOther.GetInputGeometry(), rOther.GetNodalDistances()))`

# Calls

- [GetNodalDistances](../../../functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetNodalDistances.md)