---
type: C++ Method
title: DeterminantOfJacobian
resource: geometries/quadrature_point_surface_in_volume_geometry.h#L238-L251
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/math_utils/GeneralizedDet
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`double DeterminantOfJacobian( IndexType IntegrationPointIndex, GeometryData::IntegrationMethod ThisMethod) const override`

# Calls

- [prod](../../../../functions/includes/amatrix_interface/prod.md)
- [GeneralizedDet](../../../../functions/utilities/math_utils/GeneralizedDet.md)