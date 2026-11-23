---
type: C++ Method
title: CheckPassingThroughByExtraRays
resource: modeler/internals/cartesian_ray.h#L266-L293
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/modeler/internals/cartesian_ray/CartesianRay/AddIntersection
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/internals/cartesian_ray/CartesianRay/CollapseIntersectionPoints
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool CheckPassingThroughByExtraRays(typename std::vector<std::pair<double, const TGeometryType*>>::iterator Begin, typename std::vector<std::pair<double, const TGeometryType*>>::iterator End, double Tolerance, double Delta)`

# Calls

- [AddIntersection](../../../../../functions/modeler/internals/cartesian_ray/CartesianRay/AddIntersection.md)

# Called by

- [CollapseIntersectionPoints](../../../../../functions/modeler/internals/cartesian_ray/CartesianRay/CollapseIntersectionPoints.md)