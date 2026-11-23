---
type: C++ Method
title: CreateIntegrationPoints1DGrid
resource: integration/integration_point_utilities.cpp#L118-L149
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/integration/integration_point/IntegrationPoint/Weight
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/integration/integration_point_utilities/IntegrationPointUtilities/CreateIntegrationPoints1D
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void IntegrationPointUtilities::CreateIntegrationPoints1DGrid( IntegrationPointsArrayType& rIntegrationPoints, const std::vector<double>& rSpansLocalSpace, const SizeType IntegrationPointsPerSpan)`

# Calls

- [Weight](../../../../functions/integration/integration_point/IntegrationPoint/Weight.md)

# Called by

- [CreateIntegrationPoints1D](../../../../functions/integration/integration_point_utilities/IntegrationPointUtilities/CreateIntegrationPoints1D.md)