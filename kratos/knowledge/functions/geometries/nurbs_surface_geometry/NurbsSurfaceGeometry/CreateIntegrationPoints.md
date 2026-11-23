---
type: C++ Method
title: CreateIntegrationPoints
resource: geometries/nurbs_surface_geometry.h#L557-L566
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/integration/integration_point_utilities/IntegrationPointUtilities/IntegrationPoints2D
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

- [IntegrationPoints2D](../../../../functions/integration/integration_point_utilities/IntegrationPointUtilities/IntegrationPoints2D.md)
- [GetT0](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_interval/NurbsInterval/GetT0.md)
- [GetT1](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_interval/NurbsInterval/GetT1.md)