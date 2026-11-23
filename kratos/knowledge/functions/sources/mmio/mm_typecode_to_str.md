---
type: C++ Function
title: mm_typecode_to_str
resource: sources/mmio.c#L456-L512
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/mmio/mm_strdup
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/matrix_market_interface/ReadMatrixMarketMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/matrix_market_interface/ReadMatrixMarketVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/mmio/mm_read_unsymmetric_sparse
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/mmio/mm_write_banner
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/mmio/mm_write_mtx_crd
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`char *mm_typecode_to_str(MM_typecode matcode)`

# Calls

- [mm_strdup](../../../functions/sources/mmio/mm_strdup.md)

# Called by

- [ReadMatrixMarketMatrix](../../../functions/sources/matrix_market_interface/ReadMatrixMarketMatrix.md)
- [ReadMatrixMarketVector](../../../functions/sources/matrix_market_interface/ReadMatrixMarketVector.md)
- [mm_read_unsymmetric_sparse](../../../functions/sources/mmio/mm_read_unsymmetric_sparse.md)
- [mm_write_banner](../../../functions/sources/mmio/mm_write_banner.md)
- [mm_write_mtx_crd](../../../functions/sources/mmio/mm_write_mtx_crd.md)