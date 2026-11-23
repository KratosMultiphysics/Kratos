---
type: C++ Function
title: SpMV
resource: containers/csr_matrix.h#L429-L446
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/csr_matrix/CsrMatrix/nnz
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/quaternion/y
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/quaternion/x
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SpMV(const TInputVectorType& x, TOutputVectorType& y) const`

# Calls

- [nnz](../../../functions/containers/csr_matrix/CsrMatrix/nnz.md)
- [y](../../../functions/utilities/quaternion/y.md)
- [x](../../../functions/utilities/quaternion/x.md)