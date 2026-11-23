---
type: C++ Method
title: UpdateFromLocalCoordinates
resource: utilities/quadrature_points_utility.h#L464-L489
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/GetDefaultIntegrationMethod
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/geometries/test_quadrature_point_geometry/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void UpdateFromLocalCoordinates( typename GeometryType::Pointer pGeometry, const array_1d<double, 3>& rLocalCoordinates, const double rIntegrationWeight, GeometryType& rParentGeometry)`

# Calls

- [GetDefaultIntegrationMethod](../../../../functions/geometries/geometry/GetDefaultIntegrationMethod.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/geometries/test_quadrature_point_geometry/KRATOS_TEST_CASE_IN_SUITE.md)