---
type: C++ Method
title: GetIntegrationPointData
resource: modeler/internals/skin_intersection.h#L79-L119
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/modeler/internals/skin_intersection/SkinIntersection/AddTriangleIntegrationPointData
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`double GetIntegrationPointData( std::vector<Vector>& rIntegrationPointData, const Point& rLowPoint, const Point& rHighPoint, std::size_t I, std::size_t J, std::size_t K, const GeometryData::IntegrationMethod Method, double Thickness) const`

# Calls

- [AddTriangleIntegrationPointData](../../../../../functions/modeler/internals/skin_intersection/SkinIntersection/AddTriangleIntegrationPointData.md)