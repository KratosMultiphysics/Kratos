---
type: C++ Method
title: KnotRefinementV
resource: geometries/nurbs_shape_function_utilities/nurbs_volume_refinement_utilities.h#L194-L322
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/brep_surface/BrepSurface/KnotsV
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/NumberOfControlPointsW
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_utilities/GetLowerSpan
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void KnotRefinementV(NurbsVolumeGeometryType& rGeometry, std::vector<double>& rKnotsVToInsert, PointerVector<NodeType>& rPointsRefined, Vector& rKnotsVRefined )`

# Calls

- [KnotsV](../../../../../functions/geometries/brep_surface/BrepSurface/KnotsV.md)
- [NumberOfControlPointsW](../../../../../functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/NumberOfControlPointsW.md)
- [GetLowerSpan](../../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_utilities/GetLowerSpan.md)