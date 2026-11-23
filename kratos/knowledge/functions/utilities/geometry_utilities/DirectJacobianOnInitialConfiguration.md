---
type: C++ Function
title: DirectJacobianOnInitialConfiguration
resource: utilities/geometry_utilities.h#L762-L785
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void DirectJacobianOnInitialConfiguration( GeometryType const& rGeometry, TMatrix& rJ0, const IndexType PointNumber, const GeometryType::IntegrationMethod& rIntegrationMethod )`

# Calls

- [PointsNumber](../../../functions/geometries/geometry/PointsNumber.md)