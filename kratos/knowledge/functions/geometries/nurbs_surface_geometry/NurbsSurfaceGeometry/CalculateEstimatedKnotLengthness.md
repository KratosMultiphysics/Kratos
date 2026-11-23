---
type: C++ Method
title: CalculateEstimatedKnotLengthness
resource: geometries/nurbs_surface_geometry.h#L493-L526
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_utilities/GetLowerSpan
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/brep_surface/BrepSurface/KnotsU
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/brep_surface/BrepSurface/KnotsV
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/geometries/nurbs_surface_geometry/NurbsSurfaceGeometry/Calculate
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateEstimatedKnotLengthness( CoordinatesArrayType& rKnotLengthness, const CoordinatesArrayType& rLocalCoordinates) const`

# Calls

- [GetLowerSpan](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_utilities/GetLowerSpan.md)
- [KnotsU](../../../../functions/geometries/brep_surface/BrepSurface/KnotsU.md)
- [KnotsV](../../../../functions/geometries/brep_surface/BrepSurface/KnotsV.md)
- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)

# Called by

- [Calculate](../../../../functions/geometries/nurbs_surface_geometry/NurbsSurfaceGeometry/Calculate.md)