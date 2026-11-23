---
type: C++ Method
title: KnotRefinementV
resource: geometries/nurbs_shape_function_utilities/nurbs_surface_refinement_utilities.cpp#L185-L343
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/brep_surface/BrepSurface/KnotsV
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

`void NurbsSurfaceRefinementUtilities::KnotRefinementV( NurbsSurfaceGeometryType& rGeometry, std::vector<double>& rKnotsVToInsert, PointerVector<NodeType>& rPointsRefined, Vector& rKnotsVRefined, Vector& rWeightsRefined)`

# Calls

- [KnotsV](../../../../../functions/geometries/brep_surface/BrepSurface/KnotsV.md)
- [fill](../../../../../functions/includes/amatrix_interface/Matrix/fill.md)
- [GetUpperSpan](../../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_utilities/GetUpperSpan.md)