---
type: C++ Method
title: ComputePlaneApproximation
resource: utilities/plane_approximation_utility.h#L95-L102
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/plane_approximation_utility/PlaneApproximationUtility/ComputeBasePoint
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/plane_approximation_utility/PlaneApproximationUtility/ComputePlaneNormal
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void ComputePlaneApproximation( const std::vector< array_1d<double,3> > &rPointsCoords, array_1d<double,3> &rPlaneBasePointCoords, array_1d<double,3> &rPlaneNormal)`

# Calls

- [ComputeBasePoint](../../../../functions/utilities/plane_approximation_utility/PlaneApproximationUtility/ComputeBasePoint.md)
- [ComputePlaneNormal](../../../../functions/utilities/plane_approximation_utility/PlaneApproximationUtility/ComputePlaneNormal.md)