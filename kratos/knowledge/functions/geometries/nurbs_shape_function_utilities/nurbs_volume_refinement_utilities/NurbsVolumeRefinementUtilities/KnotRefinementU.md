---
type: C++ Method
title: KnotRefinementU
resource: geometries/nurbs_shape_function_utilities/nurbs_volume_refinement_utilities.h#L54-L182
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/brep_surface/BrepSurface/KnotsU
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

`static void KnotRefinementU(NurbsVolumeGeometryType& rGeometry, std::vector<double>& rKnotsUToInsert, PointerVector<NodeType>& rPointsRefined, Vector& rKnotsURefined )`

# Calls

- [KnotsU](../../../../../functions/geometries/brep_surface/BrepSurface/KnotsU.md)
- [NumberOfControlPointsW](../../../../../functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/NumberOfControlPointsW.md)
- [GetLowerSpan](../../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_utilities/GetLowerSpan.md)