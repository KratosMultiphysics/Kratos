---
type: C++ Function
title: WriteMatrixMarketVector
resource: sources/matrix_market_interface.cpp#L571-L607
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/mmio/mm_write_banner
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/mmio/mm_write_mtx_array_size
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool WriteMatrixMarketVector(const char *FileName, const VectorType &V)`

# Calls

- [mm_write_banner](../../../functions/sources/mmio/mm_write_banner.md)
- [mm_write_mtx_array_size](../../../functions/sources/mmio/mm_write_mtx_array_size.md)