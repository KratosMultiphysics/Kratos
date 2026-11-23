---
type: C++ Method
title: ComputeQuadraticErrorFunctionPoint
resource: utilities/qef_utility.cpp#L178-L224
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/svd_utils/JacobiSingularValueDecomposition
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/qef_utility/QuadraticErrorFunction/QuadraticErrorFunctionPoint
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`array_1d<double, 3> QuadraticErrorFunction::ComputeQuadraticErrorFunctionPoint( const BoundedMatrix<double,3,3>& rATA, const BoundedMatrix<double,3,1>& rATB, const BoundedMatrix<double,3,1>& rMatCenter )`

# Calls

- [JacobiSingularValueDecomposition](../../../../functions/utilities/svd_utils/JacobiSingularValueDecomposition.md)
- [prod](../../../../functions/includes/amatrix_interface/prod.md)

# Called by

- [QuadraticErrorFunctionPoint](../../../../functions/utilities/qef_utility/QuadraticErrorFunction/QuadraticErrorFunctionPoint.md)