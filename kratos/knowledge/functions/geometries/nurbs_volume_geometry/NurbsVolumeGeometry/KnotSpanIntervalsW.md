---
type: C++ Method
title: KnotSpanIntervalsW
resource: geometries/nurbs_volume_geometry.h#L637-L654
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/NumberOfKnotsW
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/CreateIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::vector<NurbsInterval> KnotSpanIntervalsW() const`

# Calls

- [NumberOfKnotsW](../../../../functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/NumberOfKnotsW.md)

# Called by

- [CreateIntegrationPoints](../../../../functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/CreateIntegrationPoints.md)