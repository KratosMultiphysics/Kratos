---
type: C++ Method
title: OrientedBoundingBox
resource: geometries/oriented_bounding_box.cpp#L85-L150
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/math_utils/OrthonormalBasis
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/mortar_utilities/RotatePoint
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/math_utils/UnitCrossProduct
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`OrientedBoundingBox<3>::OrientedBoundingBox( const GeometryType& rGeometry, const double BoundingBoxFactor, const bool BuildFromBoundingBox )`

# Calls

- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)
- [OrthonormalBasis](../../../../functions/utilities/math_utils/OrthonormalBasis.md)
- [RotatePoint](../../../../functions/utilities/mortar_utilities/RotatePoint.md)
- [UnitCrossProduct](../../../../functions/utilities/math_utils/UnitCrossProduct.md)