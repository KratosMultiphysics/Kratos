---
type: C++ Function
title: ReadMatrixMarketVector
resource: sources/matrix_market_interface.cpp#L469-L558
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/mmio/mm_read_banner
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/mmio/mm_is_valid
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/mmio/mm_typecode_to_str
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/mmio/mm_read_mtx_array_size
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_matrix_market_interface/TestMatrixMarketInterface/test_mm_vector_io
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_matrix_market_interface/TestMatrixMarketInterface/test_mm_vector_io_cplx
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool ReadMatrixMarketVector(const char *FileName, VectorType &V)`

# Calls

- [mm_read_banner](../../../functions/sources/mmio/mm_read_banner.md)
- [mm_is_valid](../../../functions/sources/mmio/mm_is_valid.md)
- [mm_typecode_to_str](../../../functions/sources/mmio/mm_typecode_to_str.md)
- [mm_read_mtx_array_size](../../../functions/sources/mmio/mm_read_mtx_array_size.md)

# Called by

- [test_mm_vector_io](../../../functions/tests/test_matrix_market_interface/TestMatrixMarketInterface/test_mm_vector_io.md)
- [test_mm_vector_io_cplx](../../../functions/tests/test_matrix_market_interface/TestMatrixMarketInterface/test_mm_vector_io_cplx.md)