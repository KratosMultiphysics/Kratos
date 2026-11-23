---
type: C++ Function
title: ParallelProductNoAdd
resource: spaces/ublas_space.h#L1025-L1052
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/spaces/ublas_space/partial_product_no_add
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/spaces/ublas_space/Mult
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void ParallelProductNoAdd(const MatrixType& A, const VectorType& in, VectorType& out)`

# Calls

- [partial_product_no_add](../../../functions/spaces/ublas_space/partial_product_no_add.md)

# Called by

- [Mult](../../../functions/spaces/ublas_space/Mult.md)