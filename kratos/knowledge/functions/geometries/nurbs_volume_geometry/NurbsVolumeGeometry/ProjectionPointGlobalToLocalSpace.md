---
type: C++ Method
title: ProjectionPointGlobalToLocalSpace
resource: geometries/nurbs_volume_geometry.h#L986-L1009
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/brep_surface/BrepSurface/KnotsU
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/brep_surface/BrepSurface/KnotsV
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/NumberOfKnotsW
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`int ProjectionPointGlobalToLocalSpace( const CoordinatesArrayType& rPointGlobalCoordinates, CoordinatesArrayType& rProjectedPointLocalCoordinates, const double Tolerance = std::numeric_limits<double>::epsilon() ) const override`

# Calls

- [KnotsU](../../../../functions/geometries/brep_surface/BrepSurface/KnotsU.md)
- [KnotsV](../../../../functions/geometries/brep_surface/BrepSurface/KnotsV.md)
- [NumberOfKnotsW](../../../../functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/NumberOfKnotsW.md)