---
type: C++ Method
title: QuadrilateralInterface3D4
resource: geometries/quadrilateral_interface_3d_4.h#L227-L256
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

`QuadrilateralInterface3D4( typename PointType::Pointer pFirstPoint, typename PointType::Pointer pSecondPoint, typename PointType::Pointer pThirdPoint, typename PointType::Pointer pFourthPoint ) : BaseType( PointsArrayType(), &msGeometryData )`

# Calls

- [Norm3](../../../../functions/utilities/math_utils/Norm3.md)
- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)