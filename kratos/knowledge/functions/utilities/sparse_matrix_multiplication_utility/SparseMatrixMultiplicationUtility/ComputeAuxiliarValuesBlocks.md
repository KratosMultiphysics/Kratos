---
type: C++ Method
title: ComputeAuxiliarValuesBlocks
resource: utilities/sparse_matrix_multiplication_utility.h#L872-L896
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/sparse_matrix_multiplication_utility/SparseMatrixMultiplicationUtility/AssembleSparseMatrixByBlocks
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline void ComputeAuxiliarValuesBlocks( const CompressedMatrix& rMatrix, IndexType* AuxIndex2, double* AuxVals, const int CurrentRow, IndexType& RowEnd, const SizeType InitialIndexColumn, const double ContributionCoefficient = 1.0 )`

# Called by

- [AssembleSparseMatrixByBlocks](../../../../functions/utilities/sparse_matrix_multiplication_utility/SparseMatrixMultiplicationUtility/AssembleSparseMatrixByBlocks.md)