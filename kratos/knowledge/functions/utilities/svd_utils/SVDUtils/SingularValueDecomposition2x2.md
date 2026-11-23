---
type: C++ Method
title: SingularValueDecomposition2x2
resource: utilities/svd_utils.h#L312-L342
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/svd_utils/SVDUtils/SingularValueDecomposition2x2Symmetric
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/svd_utils/Jacobi
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline void SingularValueDecomposition2x2( const BoundedMatrix<TDataType, 2, 2>& rInputMatrix, BoundedMatrix<TDataType, 2, 2>& rUMatrix, BoundedMatrix<TDataType, 2, 2>& rSMatrix, BoundedMatrix<TDataType, 2, 2>& rVMatrix )`

# Calls

- [prod](../../../../functions/includes/amatrix_interface/prod.md)
- [SingularValueDecomposition2x2Symmetric](../../../../functions/utilities/svd_utils/SVDUtils/SingularValueDecomposition2x2Symmetric.md)

# Called by

- [Jacobi](../../../../functions/utilities/svd_utils/Jacobi.md)