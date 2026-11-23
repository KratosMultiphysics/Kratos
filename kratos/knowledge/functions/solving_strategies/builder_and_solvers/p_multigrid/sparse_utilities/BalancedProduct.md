---
type: C++ Function
title: BalancedProduct
resource: solving_strategies/builder_and_solvers/p_multigrid/sparse_utilities.hpp#L655-L731
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/csr_matrix/CsrMatrix/nnz
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void BalancedProduct(const typename TLHSSparse::MatrixType& rLhs, const typename TRHSSparse::VectorType& rRhs, typename TOutputSparse::VectorType& rOutput, const typename TOutputSparse::DataType Coefficient = static_cast<typename TOutputSparse::DataType>(1))`

# Calls

- [nnz](../../../../../functions/containers/csr_matrix/CsrMatrix/nnz.md)