---
type: C++ Function
title: mm_read_mtx_array_size
resource: sources/mmio.c#L221-L248
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/sources/matrix_market_interface/ReadMatrixMarketVector
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`int mm_read_mtx_array_size(FILE *f, int *M, int *N)`

# Called by

- [ReadMatrixMarketVector](../../../functions/sources/matrix_market_interface/ReadMatrixMarketVector.md)