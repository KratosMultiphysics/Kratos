---
type: C++ Function
title: ComputeVolume3DGeometryDerivative
resource: utilities/integration_utilities.h#L192-L220
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/GetDefaultIntegrationMethod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_point/IntegrationPoint/Weight
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/utilities/test_integration_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline double ComputeVolume3DGeometryDerivative( const unsigned int DerivativeNodeIndex, const unsigned int DerivativeDirectionIndex, const Geometry<TPointType>& rGeometry)`

# Calls

- [GetDefaultIntegrationMethod](../../../functions/geometries/geometry/GetDefaultIntegrationMethod.md)
- [Weight](../../../functions/integration/integration_point/IntegrationPoint/Weight.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_integration_utilities/KRATOS_TEST_CASE_IN_SUITE.md)