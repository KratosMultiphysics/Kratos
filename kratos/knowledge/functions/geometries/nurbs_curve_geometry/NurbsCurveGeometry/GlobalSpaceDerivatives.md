---
type: C++ Method
title: GlobalSpaceDerivatives
resource: geometries/nurbs_curve_geometry.h#L637-L665
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_curve_shape_functions/NurbsCurveShapeFunction/GetFirstNonzeroControlPoint
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void GlobalSpaceDerivatives( std::vector<CoordinatesArrayType>& rGlobalSpaceDerivatives, const CoordinatesArrayType& rLocalCoordinates, const SizeType DerivativeOrder) const override`

# Calls

- [GetFirstNonzeroControlPoint](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_curve_shape_functions/NurbsCurveShapeFunction/GetFirstNonzeroControlPoint.md)