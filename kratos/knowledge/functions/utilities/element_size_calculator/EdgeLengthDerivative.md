---
type: C++ Function
title: EdgeLengthDerivative
resource: utilities/element_size_calculator.h#L190-L198
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/element_size_calculator/ElementSizeCalculator<2,3>/MinimumElementSizeDerivative
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/element_size_calculator/ElementSizeCalculator<3,4>/MinimumElementSizeDerivative
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/element_size_calculator/ElementSizeCalculator<3,6>/MinimumElementSizeDerivative
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/element_size_calculator/ElementSizeCalculator<2,3>/AverageElementSizeDerivative
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/element_size_calculator/ElementSizeCalculator<3,4>/AverageElementSizeDerivative
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/element_size_calculator/ElementSizeCalculator<3,6>/AverageElementSizeDerivative
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/element_size_calculator/ElementSizeCalculator<3,8>/AverageElementSizeDerivative
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline double EdgeLengthDerivative( const unsigned int DerivativeNodeIndex, const unsigned int DerivativeDirectionIndex, const unsigned int NodeIndexA, const unsigned int NodeIndexB, const unsigned int EdgeDirection)`

# Called by

- [MinimumElementSizeDerivative](../../../functions/utilities/element_size_calculator/ElementSizeCalculator<2,3>/MinimumElementSizeDerivative.md)
- [MinimumElementSizeDerivative](../../../functions/utilities/element_size_calculator/ElementSizeCalculator<3,4>/MinimumElementSizeDerivative.md)
- [MinimumElementSizeDerivative](../../../functions/utilities/element_size_calculator/ElementSizeCalculator<3,6>/MinimumElementSizeDerivative.md)
- [AverageElementSizeDerivative](../../../functions/utilities/element_size_calculator/ElementSizeCalculator<2,3>/AverageElementSizeDerivative.md)
- [AverageElementSizeDerivative](../../../functions/utilities/element_size_calculator/ElementSizeCalculator<3,4>/AverageElementSizeDerivative.md)
- [AverageElementSizeDerivative](../../../functions/utilities/element_size_calculator/ElementSizeCalculator<3,6>/AverageElementSizeDerivative.md)
- [AverageElementSizeDerivative](../../../functions/utilities/element_size_calculator/ElementSizeCalculator<3,8>/AverageElementSizeDerivative.md)