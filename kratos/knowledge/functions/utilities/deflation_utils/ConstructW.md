---
type: C++ Function
title: ConstructW
resource: utilities/deflation_utils.h#L143-L583
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/csr_matrix/CsrMatrix/nnz
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/deflation_utils/ConstructScalarMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/deflation_utils/ExpandScalarMatrix
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void ConstructW(const std::size_t max_reduced_size, SparseMatrixType& rA, std::vector<int>& w, SparseMatrixType& deflatedA)`

# Calls

- [nnz](../../../functions/containers/csr_matrix/CsrMatrix/nnz.md)
- [ConstructScalarMatrix](../../../functions/utilities/deflation_utils/ConstructScalarMatrix.md)
- [ExpandScalarMatrix](../../../functions/utilities/deflation_utils/ExpandScalarMatrix.md)