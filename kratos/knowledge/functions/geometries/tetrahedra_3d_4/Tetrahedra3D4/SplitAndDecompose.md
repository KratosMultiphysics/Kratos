---
type: C++ Method
title: SplitAndDecompose
resource: geometries/tetrahedra_3d_4.h#L1489-L1613
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/tetrahedra_3d_4/Plane/DistanceTo
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/geometries/tetrahedra_3d_4/Tetrahedra3D4/HasIntersection
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SplitAndDecompose( const BaseType& tetra, Plane& plane, std::vector<BaseType>& inside) const`

# Calls

- [DistanceTo](../../../../functions/geometries/tetrahedra_3d_4/Plane/DistanceTo.md)

# Called by

- [HasIntersection](../../../../functions/geometries/tetrahedra_3d_4/Tetrahedra3D4/HasIntersection.md)