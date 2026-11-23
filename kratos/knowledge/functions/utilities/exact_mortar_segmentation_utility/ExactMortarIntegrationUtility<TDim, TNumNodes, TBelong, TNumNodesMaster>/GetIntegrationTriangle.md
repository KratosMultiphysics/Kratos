---
type: C++ Method
title: GetIntegrationTriangle
resource: utilities/exact_mortar_segmentation_utility.cpp#L1353-L1370
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/integration/quadrature/Quadrature/GenerateIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`Geometry<Node>::IntegrationPointsArrayType ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>::GetIntegrationTriangle()`

# Calls

- [GenerateIntegrationPoints](../../../../functions/integration/quadrature/Quadrature/GenerateIntegrationPoints.md)