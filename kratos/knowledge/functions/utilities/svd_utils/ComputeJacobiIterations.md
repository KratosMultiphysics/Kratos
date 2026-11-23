---
type: C++ Function
title: ComputeJacobiIterations
resource: utilities/svd_utils.h#L553-L626
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/svd_utils/JacobiSingularValueDecomposition
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline IndexType ComputeJacobiIterations( const TInputMatrix& rInputMatrix, TUmatrix& rUMatrix, TSMatrix& rSMatrix, TVMatrix& rVMatrix, TJ1Matrix& rJ1, TJ2Matrix& rJ2, TAuxiliaryMatrixMN& rAuxiliaryMatrixMN, TAuxiliaryMatrixM& rAuxiliaryMatrixM, TAuxiliaryMatrixN& rAuxiliaryMatrixN, const TDataType Tolerance = ZeroTolerance, const IndexType MaxIter = 200 )`

# Calls

- [prod](../../../functions/includes/amatrix_interface/prod.md)

# Called by

- [JacobiSingularValueDecomposition](../../../functions/utilities/svd_utils/JacobiSingularValueDecomposition.md)