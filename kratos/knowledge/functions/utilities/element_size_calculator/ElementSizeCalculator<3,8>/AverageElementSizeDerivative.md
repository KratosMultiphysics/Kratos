---
type: C++ Method
title: AverageElementSizeDerivative
resource: utilities/element_size_calculator.cpp#L1144-L1210
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/element_size_calculator/EdgeLengthDerivative
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`double ElementSizeCalculator<3,8>::AverageElementSizeDerivative( const unsigned int DerivativeNodeIndex, const unsigned int DerivativeDirectionIndex, const Geometry<Node >& rGeometry)`

# Calls

- [EdgeLengthDerivative](../../../../functions/utilities/element_size_calculator/EdgeLengthDerivative.md)