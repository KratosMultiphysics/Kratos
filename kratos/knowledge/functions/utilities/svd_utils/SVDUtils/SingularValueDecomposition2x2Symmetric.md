---
type: C++ Method
title: SingularValueDecomposition2x2Symmetric
resource: utilities/svd_utils.h#L357-L425
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/math_utils/Sign
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/svd_utils/SVDUtils/SingularValueDecomposition2x2
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline void SingularValueDecomposition2x2Symmetric( const BoundedMatrix<TDataType, 2, 2>& rInputMatrix, BoundedMatrix<TDataType, 2, 2>& rUMatrix, BoundedMatrix<TDataType, 2, 2>& rSMatrix, BoundedMatrix<TDataType, 2, 2>& rVMatrix )`

# Calls

- [Sign](../../../../functions/utilities/math_utils/Sign.md)
- [prod](../../../../functions/includes/amatrix_interface/prod.md)

# Called by

- [SingularValueDecomposition2x2](../../../../functions/utilities/svd_utils/SVDUtils/SingularValueDecomposition2x2.md)