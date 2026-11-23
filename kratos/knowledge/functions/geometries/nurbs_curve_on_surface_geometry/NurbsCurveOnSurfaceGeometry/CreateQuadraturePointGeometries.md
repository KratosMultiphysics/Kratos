---
type: C++ Method
title: CreateQuadraturePointGeometries
resource: geometries/nurbs_curve_on_surface_geometry.h#L592-L668
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/GetDefaultIntegrationMethod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/brep_surface/BrepSurface/KnotsU
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/brep_surface/BrepSurface/KnotsV
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CreateQuadraturePointGeometries( GeometriesArrayType& rResultGeometries, IndexType NumberOfShapeFunctionDerivatives, const IntegrationPointsArrayType& rIntegrationPoints, IntegrationInfo& rIntegrationInfo) override`

# Calls

- [GetDefaultIntegrationMethod](../../../../functions/geometries/geometry/GetDefaultIntegrationMethod.md)
- [KnotsU](../../../../functions/geometries/brep_surface/BrepSurface/KnotsU.md)
- [KnotsV](../../../../functions/geometries/brep_surface/BrepSurface/KnotsV.md)