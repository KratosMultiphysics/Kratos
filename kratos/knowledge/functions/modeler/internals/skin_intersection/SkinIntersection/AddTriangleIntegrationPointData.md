---
type: C++ Method
title: AddTriangleIntegrationPointData
resource: modeler/internals/skin_intersection.h#L245-L273
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_point/IntegrationPoint/Weight
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/internals/skin_intersection/SkinIntersection/GetIntegrationPointData
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AddTriangleIntegrationPointData( std::vector<Vector>& rIntegrationPointData, const SplitGeometryType& rGeom, const GeometryData::IntegrationMethod Method) const`

# Calls

- [PointsNumber](../../../../../functions/geometries/geometry/PointsNumber.md)
- [Weight](../../../../../functions/integration/integration_point/IntegrationPoint/Weight.md)

# Called by

- [GetIntegrationPointData](../../../../../functions/modeler/internals/skin_intersection/SkinIntersection/GetIntegrationPointData.md)