---
type: C++ Method
title: FindNode
resource: utilities/brute_force_point_locator.cpp#L23-L48
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/brute_force_point_locator/BruteForcePointLocator/NodeIsCloseEnough
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/brute_force_point_locator/BruteForcePointLocator/CheckResults
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

`int BruteForcePointLocator::FindNode(const Point& rThePoint, const Globals::Configuration configuration, const double DistanceThreshold) const`

# Calls

- [NodeIsCloseEnough](../../../../functions/utilities/brute_force_point_locator/BruteForcePointLocator/NodeIsCloseEnough.md)
- [CheckResults](../../../../functions/utilities/brute_force_point_locator/BruteForcePointLocator/CheckResults.md)

# Called by

- [__SearchPoint](../../../../functions/python_scripts/point_output_process/PointOutputProcess/__SearchPoint.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_brute_force_point_locator/KRATOS_TEST_CASE_IN_SUITE.md)