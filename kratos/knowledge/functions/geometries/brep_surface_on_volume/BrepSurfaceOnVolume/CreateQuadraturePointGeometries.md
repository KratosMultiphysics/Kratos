---
type: C++ Method
title: CreateQuadraturePointGeometries
resource: geometries/brep_surface_on_volume.h#L329-L348
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/nurbs_surface_on_volume_geometry/NurbsSurfaceOnVolumeGeometry/CreateQuadraturePointGeometriesSbm
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CreateQuadraturePointGeometries( GeometriesArrayType& rResultGeometries, IndexType NumberOfShapeFunctionDerivatives, const IntegrationPointsArrayType& rIntegrationPoints, IntegrationInfo& rIntegrationInfo) override`

# Calls

- [CreateQuadraturePointGeometriesSbm](../../../../functions/geometries/nurbs_surface_on_volume_geometry/NurbsSurfaceOnVolumeGeometry/CreateQuadraturePointGeometriesSbm.md)