---
type: C++ Method
title: ElementIsSplit
resource: utilities/embedded_skin_utility.cpp#L284-L298
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/GenerateSkin
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool inline EmbeddedSkinUtility<TDim>::ElementIsSplit( const Geometry<Node> &rGeometry, const Vector &rNodalDistances)`

# Calls

- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)

# Called by

- [GenerateSkin](../../../../functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/GenerateSkin.md)