---
type: C++ Function
title: SortRows
resource: utilities/sparse_matrix_multiplication_utility.h#L617-L647
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/sparse_matrix_multiplication_utility/MatrixMultiplicationSaad
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

`static inline void SortRows( const TIndexType* CPtr, const TSize NRows, const TSize NCols, Col* Columns, ValueType* Values )`

# Called by

- [MatrixMultiplicationSaad](../../../functions/utilities/sparse_matrix_multiplication_utility/MatrixMultiplicationSaad.md)
- [MatrixAdd](../../../functions/utilities/sparse_matrix_multiplication_utility/MatrixAdd.md)
- [TransposeMatrix](../../../functions/utilities/sparse_matrix_multiplication_utility/TransposeMatrix.md)