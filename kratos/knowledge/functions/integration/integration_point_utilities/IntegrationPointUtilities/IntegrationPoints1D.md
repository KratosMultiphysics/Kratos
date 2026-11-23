---
type: C++ Method
title: IntegrationPoints1D
resource: integration/integration_point_utilities.cpp#L151-L172
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/integration/integration_point/IntegrationPoint/Weight
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/integration/integration_point_utilities/IntegrationPointUtilities/CreateIntegrationPoints1DGauss
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void IntegrationPointUtilities::IntegrationPoints1D( typename IntegrationPointsArrayType::iterator& rIntegrationPointsBegin, SizeType PointsInU, double U0, double U1)`

# Calls

- [Weight](../../../../functions/integration/integration_point/IntegrationPoint/Weight.md)

# Called by

- [CreateIntegrationPoints1DGauss](../../../../functions/integration/integration_point_utilities/IntegrationPointUtilities/CreateIntegrationPoints1DGauss.md)