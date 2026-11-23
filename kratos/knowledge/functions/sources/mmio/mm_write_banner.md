---
type: C++ Function
title: mm_write_banner
resource: sources/mmio.c#L387-L398
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/mmio/mm_typecode_to_str
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/matrix_market_interface/WriteMatrixMarketMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/matrix_market_interface/WriteMatrixMarketVector
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`int mm_write_banner(FILE *f, MM_typecode matcode)`

# Calls

- [mm_typecode_to_str](../../../functions/sources/mmio/mm_typecode_to_str.md)

# Called by

- [WriteMatrixMarketMatrix](../../../functions/sources/matrix_market_interface/WriteMatrixMarketMatrix.md)
- [WriteMatrixMarketVector](../../../functions/sources/matrix_market_interface/WriteMatrixMarketVector.md)