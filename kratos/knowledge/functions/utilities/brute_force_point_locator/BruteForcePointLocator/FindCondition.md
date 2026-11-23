---
type: C++ Method
title: FindCondition
resource: utilities/brute_force_point_locator.cpp#L65-L78
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/brute_force_point_locator/BruteForcePointLocator/FindObject
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/point_output_process/PointOutputProcess/__SearchPoint
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_brute_force_point_locator/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`int BruteForcePointLocator::FindCondition(const Point& rThePoint, Vector& rShapeFunctionValues, const Globals::Configuration configuration, const double LocalCoordTol) const`

# Calls

- [FindObject](../../../../functions/utilities/brute_force_point_locator/BruteForcePointLocator/FindObject.md)

# Called by

- [__SearchPoint](../../../../functions/python_scripts/point_output_process/PointOutputProcess/__SearchPoint.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_brute_force_point_locator/KRATOS_TEST_CASE_IN_SUITE.md)