---
type: C++ Method
title: TriangleBoxOverlap
resource: utilities/geometry_utilities.cpp#L736-L808
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/geometry_utilities/GeometryUtils/AxisTestX
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/GeometryUtils/AxisTestY
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/GeometryUtils/PlaneBoxOverlap
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/geometries/tetrahedra_3d_4/Tetrahedra3D4/HasIntersection
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/triangle_3d_3/HasIntersection
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool GeometryUtils::TriangleBoxOverlap( const Point& rBoxCenter, const Point& rBoxHalfSize, const Point& rVertex0, const Point& rVertex1, const Point& rVertex2 )`

# Calls

- [AxisTestX](../../../../functions/utilities/geometry_utilities/GeometryUtils/AxisTestX.md)
- [AxisTestY](../../../../functions/utilities/geometry_utilities/GeometryUtils/AxisTestY.md)
- [PlaneBoxOverlap](../../../../functions/utilities/geometry_utilities/GeometryUtils/PlaneBoxOverlap.md)

# Called by

- [HasIntersection](../../../../functions/geometries/tetrahedra_3d_4/Tetrahedra3D4/HasIntersection.md)
- [HasIntersection](../../../../functions/geometries/triangle_3d_3/HasIntersection.md)