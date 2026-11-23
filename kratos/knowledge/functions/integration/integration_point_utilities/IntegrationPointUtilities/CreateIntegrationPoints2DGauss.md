---
type: C++ Method
title: CreateIntegrationPoints2DGauss
resource: integration/integration_point_utilities.cpp#L80-L110
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/integration/integration_point_utilities/IntegrationPointUtilities/IntegrationPoints2D
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/integration/integration_point_utilities/IntegrationPointUtilities/CreateIntegrationPoints2D
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void IntegrationPointUtilities::CreateIntegrationPoints2DGauss( IntegrationPointsArrayType& rIntegrationPoints, const std::vector<double>& rSpansLocalSpaceU, const std::vector<double>& rSpansLocalSpaceV, const SizeType IntegrationPointsPerSpanU, const SizeType IntegrationPointsPerSpanV)`

# Calls

- [IntegrationPoints2D](../../../../functions/integration/integration_point_utilities/IntegrationPointUtilities/IntegrationPoints2D.md)

# Called by

- [CreateIntegrationPoints2D](../../../../functions/integration/integration_point_utilities/IntegrationPointUtilities/CreateIntegrationPoints2D.md)