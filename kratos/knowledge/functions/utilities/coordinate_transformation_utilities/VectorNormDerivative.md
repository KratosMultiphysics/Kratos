---
type: C++ Function
title: VectorNormDerivative
resource: utilities/coordinate_transformation_utilities.h#L1324-L1330
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/coordinate_transformation_utilities/CalculateRotationOperatorPureShapeSensitivities
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`inline double VectorNormDerivative( const double ValueNorm, const array_1d<double, 3>& rValue, const array_1d<double, 3>& rValueDerivative) const`

# Called by

- [CalculateRotationOperatorPureShapeSensitivities](../../../functions/utilities/coordinate_transformation_utilities/CalculateRotationOperatorPureShapeSensitivities.md)