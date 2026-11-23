---
type: C++ Method
title: CreateIntegrationPoints1DGauss
resource: integration/integration_point_utilities.cpp#L59-L78
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/integration/integration_point_utilities/IntegrationPointUtilities/IntegrationPoints1D
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/integration/integration_point_utilities/IntegrationPointUtilities/CreateIntegrationPoints1D
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void IntegrationPointUtilities::CreateIntegrationPoints1DGauss( IntegrationPointsArrayType& rIntegrationPoints, const std::vector<double>& rSpansLocalSpace, const SizeType IntegrationPointsPerSpan)`

# Calls

- [IntegrationPoints1D](../../../../functions/integration/integration_point_utilities/IntegrationPointUtilities/IntegrationPoints1D.md)

# Called by

- [CreateIntegrationPoints1D](../../../../functions/integration/integration_point_utilities/IntegrationPointUtilities/CreateIntegrationPoints1D.md)