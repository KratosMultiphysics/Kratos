---
type: C++ Method
title: SpansLocalSpace
resource: geometries/nurbs_curve_on_surface_geometry.h#L246-L252
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_interval/NurbsInterval/GetT0
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_interval/NurbsInterval/GetT1
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/curve_axis_intersection/CurveAxisIntersection/ComputeAxisIntersection
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SpansLocalSpace( std::vector<double>& rSpans, IndexType DirectionIndex = 0) const override`

# Calls

- [GetT0](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_interval/NurbsInterval/GetT0.md)
- [GetT1](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_interval/NurbsInterval/GetT1.md)
- [ComputeAxisIntersection](../../../../functions/utilities/curve_axis_intersection/CurveAxisIntersection/ComputeAxisIntersection.md)