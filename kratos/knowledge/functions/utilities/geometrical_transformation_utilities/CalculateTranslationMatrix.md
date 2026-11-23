---
type: C++ Function
title: CalculateTranslationMatrix
resource: utilities/geometrical_transformation_utilities.cpp#L26-L36
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/apply_periodic_boundary_condition_process/ApplyPeriodicConditionProcess/CalculateTransformationMatrix
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateTranslationMatrix( const double Modulus, MatrixType& rMatrix, const DenseVector<double>& rDirOfTranslation )`

# Called by

- [CalculateTransformationMatrix](../../../functions/processes/apply_periodic_boundary_condition_process/ApplyPeriodicConditionProcess/CalculateTransformationMatrix.md)