---
type: C++ Function
title: Jacobi
resource: utilities/svd_utils.h#L441-L475
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/svd_utils/SVDUtils/SingularValueDecomposition2x2
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline void Jacobi( TMatrixType1& J1, TMatrixType2& J2, const TInputMatrix& rInputMatrix, const SizeType Size1, const SizeType Size2, const SizeType Index1, const SizeType Index2 )`

# Calls

- [SingularValueDecomposition2x2](../../../functions/utilities/svd_utils/SVDUtils/SingularValueDecomposition2x2.md)