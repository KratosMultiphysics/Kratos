---
type: C++ Function
title: CalculateRotationMatrix
resource: utilities/geometrical_transformation_utilities.cpp#L41-L96
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/apply_periodic_boundary_condition_process/ApplyPeriodicConditionProcess/CalculateTransformationMatrix
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateRotationMatrix( const double Theta, MatrixType& rMatrix, const DenseVector<double>& rAxisOfRotationVector, const DenseVector<double>& rCenterOfRotation )`

# Called by

- [CalculateTransformationMatrix](../../../functions/processes/apply_periodic_boundary_condition_process/ApplyPeriodicConditionProcess/CalculateTransformationMatrix.md)