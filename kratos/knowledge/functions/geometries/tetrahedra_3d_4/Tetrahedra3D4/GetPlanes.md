---
type: C++ Method
title: GetPlanes
resource: geometries/tetrahedra_3d_4.h#L1616-L1641
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/math_utils/UnitCrossProduct
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/geometries/tetrahedra_3d_4/Tetrahedra3D4/HasIntersection
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void GetPlanes(array_1d<Plane, 4>& plane) const`

# Calls

- [UnitCrossProduct](../../../../functions/utilities/math_utils/UnitCrossProduct.md)

# Called by

- [HasIntersection](../../../../functions/geometries/tetrahedra_3d_4/Tetrahedra3D4/HasIntersection.md)