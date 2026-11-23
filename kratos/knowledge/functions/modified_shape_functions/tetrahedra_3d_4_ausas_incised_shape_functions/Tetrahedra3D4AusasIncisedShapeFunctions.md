---
type: C++ Function
title: Tetrahedra3D4AusasIncisedShapeFunctions
resource: modified_shape_functions/tetrahedra_3d_4_ausas_incised_shape_functions.h#L211-L214
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetNodalDistances
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`Tetrahedra3D4AusasIncisedShapeFunctions(Tetrahedra3D4AusasIncisedShapeFunctions const& rOther) : Tetrahedra3D4AusasModifiedShapeFunctions(rOther.GetInputGeometry(), rOther.GetNodalDistances()), mExtraEdgeRatios(rOther.mExtraEdgeRatios)`

# Calls

- [GetNodalDistances](../../../functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetNodalDistances.md)