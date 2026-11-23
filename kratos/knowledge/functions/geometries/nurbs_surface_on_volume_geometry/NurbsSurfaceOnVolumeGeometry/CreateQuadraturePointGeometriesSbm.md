---
type: C++ Method
title: CreateQuadraturePointGeometriesSbm
resource: geometries/nurbs_surface_on_volume_geometry.h#L411-L523
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/GetDefaultIntegrationMethod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_utilities/GetLowerSpan
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/brep_surface/BrepSurface/KnotsU
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/brep_surface/BrepSurface/KnotsV
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/NumberOfControlPointsW
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/geometries/brep_surface_on_volume/BrepSurfaceOnVolume/CreateQuadraturePointGeometries
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CreateQuadraturePointGeometriesSbm( GeometriesArrayType& rResultGeometries, IndexType NumberOfShapeFunctionDerivatives, const IntegrationPointsArrayType& rIntegrationPoints, IntegrationInfo& rIntegrationInfo)`

# Calls

- [GetDefaultIntegrationMethod](../../../../functions/geometries/geometry/GetDefaultIntegrationMethod.md)
- [GetLowerSpan](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_utilities/GetLowerSpan.md)
- [KnotsU](../../../../functions/geometries/brep_surface/BrepSurface/KnotsU.md)
- [KnotsV](../../../../functions/geometries/brep_surface/BrepSurface/KnotsV.md)
- [NumberOfControlPointsW](../../../../functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/NumberOfControlPointsW.md)

# Called by

- [CreateQuadraturePointGeometries](../../../../functions/geometries/brep_surface_on_volume/BrepSurfaceOnVolume/CreateQuadraturePointGeometries.md)