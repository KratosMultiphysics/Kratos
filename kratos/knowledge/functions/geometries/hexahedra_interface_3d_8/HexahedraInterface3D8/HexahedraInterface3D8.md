---
type: C++ Method
title: HexahedraInterface3D8
resource: geometries/hexahedra_interface_3d_8.h#L287-L380
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/math_utils/Norm3
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`HexahedraInterface3D8( typename PointType::Pointer pPoint1, typename PointType::Pointer pPoint2, typename PointType::Pointer pPoint3, typename PointType::Pointer pPoint4, typename PointType::Pointer pPoint5, typename PointType::Pointer pPoint6, typename PointType::Pointer pPoint7, typename PointType::Pointer pPoint8 ) : BaseType( PointsArrayType(), &msGeometryData )`

# Calls

- [Norm3](../../../../functions/utilities/math_utils/Norm3.md)
- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)