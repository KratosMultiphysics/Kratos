---
type: C++ Method
title: AverageElementSizeDerivative
resource: utilities/element_size_calculator.cpp#L845-L864
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/integration_utilities/ComputeArea2DGeometryDerivative
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`double ElementSizeCalculator<2,4>::AverageElementSizeDerivative( const unsigned int DerivativeNodeIndex, const unsigned int DerivativeDirectionIndex, const Geometry<Node >& rGeometry)`

# Calls

- [ComputeArea2DGeometryDerivative](../../../../functions/utilities/integration_utilities/ComputeArea2DGeometryDerivative.md)