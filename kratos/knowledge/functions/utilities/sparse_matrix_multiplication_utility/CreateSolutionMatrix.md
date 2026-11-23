---
type: C++ Function
title: CreateSolutionMatrix
resource: utilities/sparse_matrix_multiplication_utility.h#L573-L606
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/sparse_matrix_multiplication_utility/MatrixMultiplicationSaad
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/sparse_matrix_multiplication_utility/MatrixMultiplicationRMerge
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/sparse_matrix_multiplication_utility/MatrixAdd
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/sparse_matrix_multiplication_utility/TransposeMatrix
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline void CreateSolutionMatrix( CMatrix& C, const TSize NRows, const TSize NCols, const Ptr* CPtr, const IndexType* AuxIndex2C, const ValueType* AuxValC )`

# Called by

- [MatrixMultiplicationSaad](../../../functions/utilities/sparse_matrix_multiplication_utility/MatrixMultiplicationSaad.md)
- [MatrixMultiplicationRMerge](../../../functions/utilities/sparse_matrix_multiplication_utility/MatrixMultiplicationRMerge.md)
- [MatrixAdd](../../../functions/utilities/sparse_matrix_multiplication_utility/MatrixAdd.md)
- [TransposeMatrix](../../../functions/utilities/sparse_matrix_multiplication_utility/TransposeMatrix.md)