---
type: C++ Function
title: AusasModifiedShapeFunctions
resource: modified_shape_functions/ausas_modified_shape_functions.h#L215-L217
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetNodalDistances
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`AusasModifiedShapeFunctions(AusasModifiedShapeFunctions const& rOther) : ModifiedShapeFunctions(rOther.GetInputGeometry(), rOther.GetNodalDistances())`

# Calls

- [GetNodalDistances](../../../functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetNodalDistances.md)