---
type: C++ Method
title: GlobalSpaceDerivatives
resource: geometries/nurbs_volume_geometry.h#L1058-L1103
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_volume_shape_functions/NurbsVolumeShapeFunction/GetFirstNonzeroControlPointW
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/NumberOfControlPointsW
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void GlobalSpaceDerivatives( std::vector<CoordinatesArrayType>& rGlobalSpaceDerivatives, const CoordinatesArrayType& rLocalCoordinates, const SizeType DerivativeOrder) const override`

# Calls

- [GetFirstNonzeroControlPointW](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_volume_shape_functions/NurbsVolumeShapeFunction/GetFirstNonzeroControlPointW.md)
- [NumberOfControlPointsW](../../../../functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/NumberOfControlPointsW.md)