---
type: C++ Function
title: ExpandScalarMatrix
resource: utilities/deflation_utils.h#L526-L556
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/deflation_utils/ConstructW
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void ExpandScalarMatrix(const std::size_t n_row, const std::size_t block_size, const std::size_t Ap[], const std::size_t Aj[], SparseMatrixType& Aexpanded )`

# Called by

- [ConstructW](../../../functions/utilities/deflation_utils/ConstructW.md)