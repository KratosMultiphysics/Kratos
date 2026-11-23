---
type: C++ Method
title: ComputeAxisIntersectionSBM
resource: geometries/nurbs_curve_on_surface_geometry.h#L331-L421
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/geometries/nurbs_curve_on_surface_geometry/NurbsCurveOnSurfaceGeometry/SpansLocalSpaceSBM
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ComputeAxisIntersectionSBM( std::vector<double>& rSpans, double Start, double End, const std::vector<double>& surface_spans_u, const std::vector<double>& surface_spans_v) const`

# Calls

- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)

# Called by

- [SpansLocalSpaceSBM](../../../../functions/geometries/nurbs_curve_on_surface_geometry/NurbsCurveOnSurfaceGeometry/SpansLocalSpaceSBM.md)