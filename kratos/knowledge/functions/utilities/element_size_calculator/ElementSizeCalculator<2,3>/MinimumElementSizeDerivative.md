---
type: C++ Method
title: MinimumElementSizeDerivative
resource: utilities/element_size_calculator.cpp#L63-L125
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

`double ElementSizeCalculator<2,3>::MinimumElementSizeDerivative( const unsigned int DerivativeNodeIndex, const unsigned int DerivativeDirectionIndex, const Geometry<Node >& rGeometry)`

# Calls

- [EdgeLengthDerivative](../../../../functions/utilities/element_size_calculator/EdgeLengthDerivative.md)