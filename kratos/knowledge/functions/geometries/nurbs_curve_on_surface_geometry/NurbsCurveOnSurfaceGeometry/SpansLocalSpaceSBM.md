---
type: C++ Method
title: SpansLocalSpaceSBM
resource: geometries/nurbs_curve_on_surface_geometry.h#L259-L265
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
  - target: functions/geometries/nurbs_curve_on_surface_geometry/NurbsCurveOnSurfaceGeometry/ComputeAxisIntersectionSBM
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SpansLocalSpaceSBM( std::vector<double>& rSpans, IndexType DirectionIndex = 0) const`

# Calls

- [GetT0](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_interval/NurbsInterval/GetT0.md)
- [GetT1](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_interval/NurbsInterval/GetT1.md)
- [ComputeAxisIntersectionSBM](../../../../functions/geometries/nurbs_curve_on_surface_geometry/NurbsCurveOnSurfaceGeometry/ComputeAxisIntersectionSBM.md)