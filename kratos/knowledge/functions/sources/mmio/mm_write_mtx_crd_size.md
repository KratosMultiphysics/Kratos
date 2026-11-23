---
type: C++ Function
title: mm_write_mtx_crd_size
resource: sources/mmio.c#L182-L188
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/sources/matrix_market_interface/WriteMatrixMarketMatrix
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`int mm_write_mtx_crd_size(FILE *f, int M, int N, int nz)`

# Called by

- [WriteMatrixMarketMatrix](../../../functions/sources/matrix_market_interface/WriteMatrixMarketMatrix.md)