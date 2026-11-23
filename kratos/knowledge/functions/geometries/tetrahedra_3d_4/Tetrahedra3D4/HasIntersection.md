---
type: C++ Method
title: HasIntersection
resource: geometries/tetrahedra_3d_4.h#L1415-L1444
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/tetrahedra_3d_4/Tetrahedra3D4/GetPlanes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/tetrahedra_3d_4/Tetrahedra3D4/SplitAndDecompose
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/GeometryUtils/TriangleBoxOverlap
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool HasIntersection( const BaseType& rThisGeometry) const override`

# Calls

- [GetPlanes](../../../../functions/geometries/tetrahedra_3d_4/Tetrahedra3D4/GetPlanes.md)
- [SplitAndDecompose](../../../../functions/geometries/tetrahedra_3d_4/Tetrahedra3D4/SplitAndDecompose.md)
- [TriangleBoxOverlap](../../../../functions/utilities/geometry_utilities/GeometryUtils/TriangleBoxOverlap.md)