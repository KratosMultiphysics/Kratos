---
type: C++ Method
title: SetNumberOfIntegrationPointsPerSpan
resource: integration/integration_info.cpp#L109-L113
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/tests/cpp_tests/geometries/test_coupling_geometry_quadrature/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_nurbs_curve/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void IntegrationInfo::SetNumberOfIntegrationPointsPerSpan(IndexType DimensionIndex, SizeType NumberOfIntegrationPointsPerSpan)`

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/geometries/test_coupling_geometry_quadrature/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/geometries/test_nurbs_curve/KRATOS_TEST_CASE_IN_SUITE.md)