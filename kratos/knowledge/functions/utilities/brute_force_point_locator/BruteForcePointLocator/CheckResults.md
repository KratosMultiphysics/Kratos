---
type: C++ Method
title: CheckResults
resource: utilities/brute_force_point_locator.cpp#L137-L150
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/brute_force_point_locator/BruteForcePointLocator/FindNode
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/brute_force_point_locator/BruteForcePointLocator/FindObject
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void BruteForcePointLocator::CheckResults(const std::string& rObjectName, const Point& rThePoint, const int LocalObjectFound) const`

# Called by

- [FindNode](../../../../functions/utilities/brute_force_point_locator/BruteForcePointLocator/FindNode.md)
- [FindObject](../../../../functions/utilities/brute_force_point_locator/BruteForcePointLocator/FindObject.md)