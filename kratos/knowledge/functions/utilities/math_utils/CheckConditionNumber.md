---
type: C++ Function
title: CheckConditionNumber
resource: utilities/math_utils.h#L252-L278
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/norm_frobenius
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/includes/mortar_classes/DualLagrangeMultiplierOperators/CalculateAe
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/mortar_classes/DualLagrangeMultiplierOperatorsWithDerivatives/CalculateDeltaAe
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/math_utils/InvertMatrix
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline bool CheckConditionNumber( const TMatrix1& rInputMatrix, TMatrix2& rInvertedMatrix, const double Tolerance = ZeroTolerance, const bool ThrowError = true )`

# Calls

- [norm_frobenius](../../../functions/includes/amatrix_interface/norm_frobenius.md)

# Called by

- [CalculateAe](../../../functions/includes/mortar_classes/DualLagrangeMultiplierOperators/CalculateAe.md)
- [CalculateDeltaAe](../../../functions/includes/mortar_classes/DualLagrangeMultiplierOperatorsWithDerivatives/CalculateDeltaAe.md)
- [InvertMatrix](../../../functions/utilities/math_utils/InvertMatrix.md)