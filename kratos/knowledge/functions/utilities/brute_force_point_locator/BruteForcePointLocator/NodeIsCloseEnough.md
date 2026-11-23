---
type: C++ Method
title: NodeIsCloseEnough
resource: utilities/brute_force_point_locator.cpp#L152-L168
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/brute_force_point_locator/BruteForcePointLocator/FindNode
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool BruteForcePointLocator::NodeIsCloseEnough(const Node& rNode, const Point& rThePoint, const Globals::Configuration configuration, const double DistanceThreshold) const`

# Called by

- [FindNode](../../../../functions/utilities/brute_force_point_locator/BruteForcePointLocator/FindNode.md)