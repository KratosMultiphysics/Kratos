---
type: C++ Method
title: FindObject
resource: utilities/brute_force_point_locator.cpp#L81-L135
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/brute_force_point_locator/BruteForcePointLocator/CheckResults
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/brute_force_point_locator/BruteForcePointLocator/FindElement
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/brute_force_point_locator/BruteForcePointLocator/FindCondition
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void BruteForcePointLocator::FindObject(const TObjectType& rObjects, const std::string& rObjectName, const Point& rThePoint, int& rObjectId, Vector& rShapeFunctionValues, const Globals::Configuration configuration, const double LocalCoordTol) const`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [CheckResults](../../../../functions/utilities/brute_force_point_locator/BruteForcePointLocator/CheckResults.md)

# Called by

- [FindElement](../../../../functions/utilities/brute_force_point_locator/BruteForcePointLocator/FindElement.md)
- [FindCondition](../../../../functions/utilities/brute_force_point_locator/BruteForcePointLocator/FindCondition.md)