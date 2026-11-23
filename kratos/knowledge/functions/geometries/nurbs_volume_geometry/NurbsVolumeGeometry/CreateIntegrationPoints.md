---
type: C++ Method
title: CreateIntegrationPoints
resource: geometries/nurbs_volume_geometry.h#L676-L686
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/KnotSpanIntervalsW
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_point_utilities/IntegrationPointUtilities/IntegrationPoints3D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_interval/NurbsInterval/GetT0
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_interval/NurbsInterval/GetT1
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CreateIntegrationPoints( IntegrationPointsArrayType& rIntegrationPoints, IntegrationInfo& rIntegrationInfo) const override`

# Calls

- [KnotSpanIntervalsW](../../../../functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/KnotSpanIntervalsW.md)
- [IntegrationPoints3D](../../../../functions/integration/integration_point_utilities/IntegrationPointUtilities/IntegrationPoints3D.md)
- [GetT0](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_interval/NurbsInterval/GetT0.md)
- [GetT1](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_interval/NurbsInterval/GetT1.md)