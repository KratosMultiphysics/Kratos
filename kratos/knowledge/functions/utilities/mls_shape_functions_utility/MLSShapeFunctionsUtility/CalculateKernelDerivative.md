---
type: C++ Method
title: CalculateKernelDerivative
resource: utilities/mls_shape_functions_utility.cpp#L43-L64
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MLSShapeFunctionsUtility::CalculateKernelDerivative( const array_1d<double,3>& rRadVect, const double h, array_1d<double,TDim>& rKernelDerivative)`

# Calls

- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)