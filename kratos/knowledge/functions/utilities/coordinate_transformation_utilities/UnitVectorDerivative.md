---
type: C++ Function
title: UnitVectorDerivative
resource: utilities/coordinate_transformation_utilities.h#L1332-L1340
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/coordinate_transformation_utilities/CalculateRotationOperatorPureShapeSensitivities
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`inline array_1d<double, 3> UnitVectorDerivative( const double VectorNorm, const double VectorNormDerivative, const array_1d<double, 3>& rVector, const array_1d<double, 3>& rVectorDerivative) const`

# Called by

- [CalculateRotationOperatorPureShapeSensitivities](../../../functions/utilities/coordinate_transformation_utilities/CalculateRotationOperatorPureShapeSensitivities.md)