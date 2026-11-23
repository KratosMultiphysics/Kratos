---
type: C++ Method
title: CreateFromLocalCoordinates
resource: utilities/quadrature_points_utility.h#L69-L102
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/GetDefaultIntegrationMethod
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/quadrature_points_utility/CreateQuadraturePointsUtility/CreateFromCoordinates
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static typename GeometryType::Pointer CreateFromLocalCoordinates( GeometryType& rGeometry, const array_1d<double, 3>& rLocalCoordinates, double integration_weight)`

# Calls

- [GetDefaultIntegrationMethod](../../../../functions/geometries/geometry/GetDefaultIntegrationMethod.md)

# Called by

- [CreateFromCoordinates](../../../../functions/utilities/quadrature_points_utility/CreateQuadraturePointsUtility/CreateFromCoordinates.md)