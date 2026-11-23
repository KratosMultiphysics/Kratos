---
type: C++ Method
title: CreateQuadraturePointGeometries
resource: geometries/nurbs_volume_geometry.h#L787-L805
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/integration/integration_info/GetQuadratureMethod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/GetDefaultIntegrationMethod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/NumberOfControlPointsW
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CreateQuadraturePointGeometries( GeometriesArrayType& rResultGeometries, IndexType NumberOfShapeFunctionDerivatives, IntegrationInfo& rIntegrationInfo) override`

# Calls

- [GetQuadratureMethod](../../../../functions/integration/integration_info/GetQuadratureMethod.md)
- [GetDefaultIntegrationMethod](../../../../functions/geometries/geometry/GetDefaultIntegrationMethod.md)
- [NumberOfControlPointsW](../../../../functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/NumberOfControlPointsW.md)