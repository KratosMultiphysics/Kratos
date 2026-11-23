---
type: C++ Function
title: Det2
resource: utilities/math_utils.h#L565-L570
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/element_size_calculator/ElementSizeCalculator<2,9>/AverageElementSizeDerivative
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/integration_utilities/ComputeArea2DGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline double Det2(const TMatrixType& rA)`

# Called by

- [AverageElementSizeDerivative](../../../functions/utilities/element_size_calculator/ElementSizeCalculator<2,9>/AverageElementSizeDerivative.md)
- [ComputeArea2DGeometry](../../../functions/utilities/integration_utilities/ComputeArea2DGeometry.md)