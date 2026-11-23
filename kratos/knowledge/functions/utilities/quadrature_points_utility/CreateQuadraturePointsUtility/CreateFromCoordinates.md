---
type: C++ Method
title: CreateFromCoordinates
resource: utilities/quadrature_points_utility.h#L54-L67
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/quadrature_points_utility/CreateQuadraturePointsUtility/CreateFromLocalCoordinates
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static typename GeometryType::Pointer CreateFromCoordinates( typename GeometryType::Pointer pGeometry, const array_1d<double, 3>& rCoordinates, double integration_weight)`

# Calls

- [CreateFromLocalCoordinates](../../../../functions/utilities/quadrature_points_utility/CreateQuadraturePointsUtility/CreateFromLocalCoordinates.md)