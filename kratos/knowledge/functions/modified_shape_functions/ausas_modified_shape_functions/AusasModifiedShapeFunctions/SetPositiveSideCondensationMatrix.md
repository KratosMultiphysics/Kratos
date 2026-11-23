---
type: C++ Method
title: SetPositiveSideCondensationMatrix
resource: modified_shape_functions/ausas_modified_shape_functions.cpp#L58-L93
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetNodalDistances
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AusasModifiedShapeFunctions::SetPositiveSideCondensationMatrix( Matrix& rPosSideCondMatrix, const std::vector<int>& rEdgeNodeI, const std::vector<int>& rEdgeNodeJ, const std::vector<int>& rSplitEdges)`

# Calls

- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)
- [GetNodalDistances](../../../../functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetNodalDistances.md)