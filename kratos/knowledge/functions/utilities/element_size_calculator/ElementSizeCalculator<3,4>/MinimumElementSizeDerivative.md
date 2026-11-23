---
type: C++ Method
title: MinimumElementSizeDerivative
resource: utilities/element_size_calculator.cpp#L347-L446
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/element_size_calculator/EdgeLengthDerivative
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/element_size_calculator/HsqDerivative3D
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`double ElementSizeCalculator<3,4>::MinimumElementSizeDerivative( const unsigned int DerivativeNodeIndex, const unsigned int DerivativeDirectionIndex, const Geometry<Node >& rGeometry)`

# Calls

- [EdgeLengthDerivative](../../../../functions/utilities/element_size_calculator/EdgeLengthDerivative.md)
- [HsqDerivative3D](../../../../functions/utilities/element_size_calculator/HsqDerivative3D.md)