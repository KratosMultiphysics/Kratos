---
type: C++ Method
title: CreateQuadraturePointGeometries
resource: geometries/coupling_geometry.h#L433-L524
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CreateQuadraturePointGeometries( GeometriesArrayType& rResultGeometries, IndexType NumberOfShapeFunctionDerivatives, const IntegrationPointsArrayType& rIntegrationPoints, IntegrationInfo& rIntegrationInfo) override`

# Calls

- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)