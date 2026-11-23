---
type: C++ Method
title: SetCondensationMatrix
resource: modified_shape_functions/modified_shape_functions.cpp#L449-L484
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModifiedShapeFunctions::SetCondensationMatrix( Matrix& rIntPointCondMatrix, const std::vector<int>& rEdgeNodeI, const std::vector<int>& rEdgeNodeJ, const std::vector<int>& rSplitEdges)`

# Calls

- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)