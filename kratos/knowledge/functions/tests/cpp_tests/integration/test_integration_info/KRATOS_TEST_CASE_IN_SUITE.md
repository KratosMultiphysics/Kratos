---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/integration/test_integration_info.cpp#L19-L30
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
  - target: functions/integration/integration_info/IntegrationInfo/SetQuadratureMethod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/flags/IsDefined
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(CheckIntegrationInfo, KratosCoreFastSuite)`

# Calls

- [GetNumberOfIntegrationPointsPerSpan](../../../../../functions/integration/integration_info/IntegrationInfo/GetNumberOfIntegrationPointsPerSpan.md)
- [GetQuadratureMethod](../../../../../functions/integration/integration_info/GetQuadratureMethod.md)
- [SetQuadratureMethod](../../../../../functions/integration/integration_info/IntegrationInfo/SetQuadratureMethod.md)
- [IsDefined](../../../../../functions/containers/flags/IsDefined.md)