---
type: C++ Function
title: Tetrahedra3D4AusasModifiedShapeFunctions
resource: modified_shape_functions/tetrahedra_3d_4_ausas_modified_shape_functions.h#L182-L189
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetNodalDistances
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`Tetrahedra3D4AusasModifiedShapeFunctions(Tetrahedra3D4AusasModifiedShapeFunctions const& rOther) : AusasModifiedShapeFunctions(rOther.GetInputGeometry(), rOther.GetNodalDistances()), mpTetrahedraSplitter(new DivideTetrahedra3D4<Node>(*rOther.GetInputGeometry(), rOther.GetNodalDistances()))`

# Calls

- [GetNodalDistances](../../../functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetNodalDistances.md)