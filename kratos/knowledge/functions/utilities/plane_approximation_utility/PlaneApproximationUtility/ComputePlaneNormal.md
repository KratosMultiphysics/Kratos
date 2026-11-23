---
type: C++ Method
title: ComputePlaneNormal
resource: utilities/plane_approximation_utility.h#L189-L215
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/plane_approximation_utility/PlaneApproximationUtility/SetMatrixA
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/math_utils/GaussSeidelEigenSystem
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/plane_approximation_utility/PlaneApproximationUtility/ComputePlaneApproximation
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void ComputePlaneNormal( const std::vector< array_1d<double,3> > &rPointsCoords, const array_1d<double,3> &rPlaneBasePointCoords, array_1d<double,3> &rPlaneNormal)`

# Calls

- [SetMatrixA](../../../../functions/utilities/plane_approximation_utility/PlaneApproximationUtility/SetMatrixA.md)
- [GaussSeidelEigenSystem](../../../../functions/utilities/math_utils/GaussSeidelEigenSystem.md)

# Called by

- [ComputePlaneApproximation](../../../../functions/utilities/plane_approximation_utility/PlaneApproximationUtility/ComputePlaneApproximation.md)