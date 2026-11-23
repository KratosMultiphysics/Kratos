---
type: C++ Method
title: PointDistanceToQuadrilateral3D
resource: utilities/geometry_utilities.cpp#L399-L410
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/geometries/hexahedra_3d_8/Hexahedra3D8/CalculateDistance
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/prism_3d_6/Prism3D6/CalculateDistance
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/pyramid_3d_5/Pyramid3D5/CalculateDistance
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/quadrilateral_3d_4/CalculateDistance
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`double GeometryUtils::PointDistanceToQuadrilateral3D( const Point& rQuadrilateralPoint1, const Point& rQuadrilateralPoint2, const Point& rQuadrilateralPoint3, const Point& rQuadrilateralPoint4, const Point& rPoint )`

# Called by

- [CalculateDistance](../../../../functions/geometries/hexahedra_3d_8/Hexahedra3D8/CalculateDistance.md)
- [CalculateDistance](../../../../functions/geometries/prism_3d_6/Prism3D6/CalculateDistance.md)
- [CalculateDistance](../../../../functions/geometries/pyramid_3d_5/Pyramid3D5/CalculateDistance.md)
- [CalculateDistance](../../../../functions/geometries/quadrilateral_3d_4/CalculateDistance.md)