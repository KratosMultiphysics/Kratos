---
type: C++ Method
title: SetMatrixA
resource: utilities/plane_approximation_utility.h#L162-L181
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/plane_approximation_utility/PlaneApproximationUtility/ComputePlaneNormal
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void SetMatrixA( const std::vector< array_1d< double,3 > > &rPointsCoords, const array_1d<double,3> &rPlaneBasePointCoords, BoundedMatrix<double,TDim,TDim> &rA)`

# Called by

- [ComputePlaneNormal](../../../../functions/utilities/plane_approximation_utility/PlaneApproximationUtility/ComputePlaneNormal.md)