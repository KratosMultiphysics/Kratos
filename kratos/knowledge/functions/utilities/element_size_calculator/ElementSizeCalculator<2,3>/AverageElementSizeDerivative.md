---
type: C++ Method
title: AverageElementSizeDerivative
resource: utilities/element_size_calculator.cpp#L777-L807
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

`double ElementSizeCalculator<2,3>::AverageElementSizeDerivative( const unsigned int DerivativeNodeIndex, const unsigned int DerivativeDirectionIndex, const Geometry<Node >& rGeometry)`

# Calls

- [EdgeLengthDerivative](../../../../functions/utilities/element_size_calculator/EdgeLengthDerivative.md)