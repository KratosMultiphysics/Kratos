---
type: C++ Method
title: ComputeNonZeroBlocks
resource: utilities/sparse_matrix_multiplication_utility.h#L846-L861
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/sparse_matrix_multiplication_utility/SparseMatrixMultiplicationUtility/AssembleSparseMatrixByBlocks
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline void ComputeNonZeroBlocks( const CompressedMatrix& rMatrix, const int CurrentRow, IndexType& rNonZeroColsAux2 )`

# Called by

- [AssembleSparseMatrixByBlocks](../../../../functions/utilities/sparse_matrix_multiplication_utility/SparseMatrixMultiplicationUtility/AssembleSparseMatrixByBlocks.md)