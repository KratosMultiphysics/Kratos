---
type: C++ Method
title: KnotRefinementU
resource: geometries/nurbs_shape_function_utilities/nurbs_surface_refinement_utilities.cpp#L23-L183
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/brep_surface/BrepSurface/KnotsU
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/Matrix/fill
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_utilities/GetUpperSpan
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void NurbsSurfaceRefinementUtilities::KnotRefinementU( NurbsSurfaceGeometryType& rGeometry, std::vector<double>& rKnotsUToInsert, PointerVector<NodeType>& rPointsRefined, Vector& rKnotsURefined, Vector& rWeightsRefined)`

# Calls

- [KnotsU](../../../../../functions/geometries/brep_surface/BrepSurface/KnotsU.md)
- [fill](../../../../../functions/includes/amatrix_interface/Matrix/fill.md)
- [GetUpperSpan](../../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_utilities/GetUpperSpan.md)