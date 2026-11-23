---
type: C++ Method
title: GetIntegrationMethod
resource: integration/integration_info.cpp#L124-L130
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/integration/integration_info/IntegrationInfo/GetNumberOfIntegrationPointsPerSpan
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_info/GetQuadratureMethod
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`IntegrationInfo::IntegrationMethod IntegrationInfo::GetIntegrationMethod( IndexType DimensionIndex) const`

# Calls

- [GetNumberOfIntegrationPointsPerSpan](../../../../functions/integration/integration_info/IntegrationInfo/GetNumberOfIntegrationPointsPerSpan.md)
- [GetQuadratureMethod](../../../../functions/integration/integration_info/GetQuadratureMethod.md)