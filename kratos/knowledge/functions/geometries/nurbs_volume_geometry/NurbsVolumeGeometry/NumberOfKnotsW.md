---
type: C++ Method
title: NumberOfKnotsW
resource: geometries/nurbs_volume_geometry.h#L360-L363
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_volume_refinement_utilities/NurbsVolumeRefinementUtilities/KnotRefinementW
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/NumberOfControlPointsW
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/DomainIntervalW
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/KnotSpanIntervalsW
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/ProjectionPointGlobalToLocalSpace
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`SizeType NumberOfKnotsW() const`

# Called by

- [KnotRefinementW](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_volume_refinement_utilities/NurbsVolumeRefinementUtilities/KnotRefinementW.md)
- [NumberOfControlPointsW](../../../../functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/NumberOfControlPointsW.md)
- [DomainIntervalW](../../../../functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/DomainIntervalW.md)
- [KnotSpanIntervalsW](../../../../functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/KnotSpanIntervalsW.md)
- [ProjectionPointGlobalToLocalSpace](../../../../functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/ProjectionPointGlobalToLocalSpace.md)